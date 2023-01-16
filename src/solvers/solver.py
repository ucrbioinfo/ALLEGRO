import sys
import numpy
import itertools
from timeit import default_timer
from ortools.linear_solver import pywraplp

from cover_set_parsers.coversets import Coversets


# For some reason, Google ORTools does not like its objects being passed around.
# So, everything goes in this class.
class Solver:
    def __init__(
        self, 
        coverset_parser: Coversets, 
        beta: int, 
        solver_engine: str, 
        exhaustive_threshold: int,
        num_trials: int,
        objective: str='max') -> None:
        
        self.coversets = coverset_parser.get_coversets()  # TODO these should be a class 
        self.species = coverset_parser.get_species_set()
        self.k = len(self.coversets)
        self.exhaustive_threshold = exhaustive_threshold
        self.solved_with_exhaustive: bool = False
        self.num_trials = num_trials
        self.objective = objective
        
        self.solver = pywraplp.Solver.CreateSolver(solver_engine)

        self.beta: int = beta

        self.set_size_for_each_trial: list[int] = list()
        self.average_score_for_each_trial: list[float] = list()
        self.num_while_iters_for_each_trial: list[int] = list()
        self.fractional_glop_vals = list()

        self.winner_scores: list[float] = list()
        self.num_feasible_guides_with_prob_lt_one: int = 0
        self.num_non_zero_feasible: int = 0
        self.num_exhausted_combos: int = 0

        self.solver_time: float = 0.0
        self.exhaustive_time: float = 0.0
        self.randomized_rounding_time: float = 0.0

        
    def solve(self) -> set[str]:
        print('Solver working...')
        
        cover: dict[int, list[str]] = dict()  # { species_id: [list of guide seqs that cover it] }
        max_score: float = numpy.NINF  # Used to transform the efficiency scores.

        for guide_seq, tupe in self.coversets.items():
            for species in tupe[1]:
                cover.setdefault(species, list()).append(guide_seq)

                if tupe[0] > max_score:
                    max_score = tupe[0]
        # -------------------------------------------


        # ------ VARIABLES --------------------------
        print('Creating variables...')

        vars = dict((seq, self.solver.NumVar(lb=0.0, ub=1.0, name=seq)) for seq in self.coversets.keys())

        print('Number of variables (guides) =', self.solver.NumVariables())
        # -------------------------------------------


        # ------ CONSTRAINTS ------------------------
        print('Creating constraints...')

        # All species must be covered by at least 1 guide
        for species in self.species:
            self.solver.Add(self.solver.Sum([vars[seq] for seq in cover[species]]) >= 1)

        print('Number of constraints (species) =', self.solver.NumConstraints())
        # -------------------------------------------


        # ------ OBJECTIVE --------------------------
        print('Creating objectives...')

        objective_1_terms = list()

        if self.objective == 'max':
            objective_2_terms = list()

            for seq, var in vars.items():
                objective_1_terms.append(1 * var)
                objective_2_terms.append(self.coversets[seq][0] * var)

            objective_1 = self.solver.Sum(objective_1_terms)
            objective_2 = self.solver.Sum(objective_2_terms)

            self.solver.Add(objective_1 <= self.beta)

            self.solver.Maximize(objective_2)

        elif self.objective == 'min':
            print('Minimizing set size (unweighted set cover). Ignoring beta...')

            for seq, var in vars.items():
                objective_1_terms.append(1 * var)
            
            objective_1 = self.solver.Sum(objective_1_terms)

            self.solver.Minimize(objective_1)

        # -------- SOLVE ----------------------------
        print('Solving...')

        solutions = list()

        start_time = default_timer()

        status = self.solver.Solve()

        self.solver_time = default_timer() - start_time
        print('\nSolver took {t} seconds.'.format(t=round(self.solver_time, 4)))

        # self.solver_wall_time = self.solver.wall_time() # solver wall time
        # print(('Problem solved in {t:.2f} milliseconds'.format(t=self.solver_wall_time)))
        # print(self.solver.ExportModelAsLpFormat(False)) # Logs variable bounds, constraints, and obj terms

        match status:
            case pywraplp.Solver.OPTIMAL | pywraplp.Solver.FEASIBLE:
                if status == pywraplp.Solver.OPTIMAL: print('OPTIMAL')
                if status == pywraplp.Solver.FEASIBLE: print('FEASIBLE')

                print('Total cost = {cost}'.format(cost=self.solver.Objective().Value()))
                print('Feasible guides (variables with a probability higher than 0.0):\n')

                for seq, var in vars.items():
                    if var.solution_value() > 0.0:
                        print(
                            'Guide: {seq} \t with probability: {p:.2f} \t with score: {score} \t covering species: {species}'.format(
                                seq=seq,
                                p=var.solution_value(),
                                score=self.coversets[seq][0],
                                species=self.coversets[seq][1],
                            ).expandtabs(10)
                        )
                        self.fractional_glop_vals.append(var.solution_value()) 
                        solutions.append(var)
                self.num_non_zero_feasible = len(solutions)
            
            case pywraplp.Solver.INFEASIBLE:
                print('Problem is proven to be infeasible.')
                return set()
            
            case pywraplp.Solver.UNBOUNDED:
                print('Problem is proven to be unbounded.')
                return set()
            
            case pywraplp.Solver.ABNORMAL:
                print('Abnormal status. Something went wrong. Call Amir.')
                sys.exit(1)

            case _:
                print('You should not have gotten here. Something went really wrong in solver.py Call Amir.')
                sys.exit(1)

        print()
        # -------------------------------------------

        for var in solutions:
            if var.solution_value() < 1:
                self.num_feasible_guides_with_prob_lt_one += 1


        winners: set[str] = set()
        len_winners: int = len(self.species)  # Worst case, have to pick a guide for every species

        len_solutions = len(solutions)
        print('Number of feasible candidate guides: {n}'.format(n=len_solutions))

        # -------- RANDOMIZED ROUND ----------------
        if len_solutions > self.exhaustive_threshold:
            print('There are more feasible guides ({n}) than the chosen threshold ({t}):'.format(
                n=len_solutions,
                t=self.exhaustive_threshold,
                )
            )
            print('Using randomized rounding with {n} trials...'.format(n=self.num_trials))
            
            rng = numpy.random.default_rng()
            trial_with_smallest_size: int = 0

            # Heuristic?
            # solutions.sort(key=lambda var: var.solution_value(), reverse=True)
            
            start_time = default_timer()

            for trial in range(1, self.num_trials + 1):
                I_this_trial: set[int] = set()
                winners_this_trial: set[str] = set()
                iters_this_trial: int = 0

                while I_this_trial != self.species:
                    iters_this_trial += 1
                    score_sum_this_trial: float = 0.0

                    for var in solutions:
                        if var.solution_value() == 1.0 or var.solution_value() > rng.random():
                            if self.coversets[var.name()][1] - I_this_trial:
                                I_this_trial = I_this_trial.union(self.coversets[var.name()][1])
                                winners_this_trial = winners_this_trial.union({var.name()})
                # print('Covered all species')

                score_sum_this_trial = 0
                for w in winners_this_trial:
                    score_sum_this_trial += self.coversets[w][0]

                len_winners_this_trial = len(winners_this_trial)
                self.set_size_for_each_trial.append(len_winners_this_trial)
                self.average_score_for_each_trial.append(score_sum_this_trial / len_winners_this_trial)
                self.num_while_iters_for_each_trial.append(iters_this_trial)
                # ------------------------------------------------------------

                # Keep the smallest set of winners
                if len_winners_this_trial <= len_winners:
                    winners = winners_this_trial
                    len_winners = len_winners_this_trial
                    trial_with_smallest_size = trial

                # Uncomment to see trial count and best size so far.
                # if trial % 100 == 0:
                #     print('Trial {t} -- Smallest solution set size so far: {b}'.format(t=trial, b=len_winners))


            self.randomized_rounding_time = default_timer() - start_time
            print('\nApproximate rounding took {t} seconds for {n} trials.'.format(t=round(self.randomized_rounding_time, 4), n=self.num_trials))
            print('This run\'s smallest solution size found last in trial {t}'.format(t=trial_with_smallest_size))
            # --------------------------------------------------------------------------------------
        
        # Exhaustive search
        else:
            self.solved_with_exhaustive = True
            print('There are fewer feasible solutions ({f}) than the chosen threshold ({t}):'.format(
                f=len_solutions, 
                t=self.exhaustive_threshold,
                )
            )
            print('Using exhaustive search...\n')
            start_time = default_timer()

            possibilities_exhausted = 0

            # ---- Keep guaranteed guides (where prob == 1) - Do not create combinations with these --------
            small_solutions = list()
            species_guaranteed_to_be_covered: set[int] = set()
            guides_guaranteed_to_be_winners: set[str] = set()

            for var in solutions:
                if var.solution_value() == 1.0:
                    if self.coversets[var.name()][1] - species_guaranteed_to_be_covered:
                        species_guaranteed_to_be_covered = species_guaranteed_to_be_covered.union(self.coversets[var.name()][1])
                        guides_guaranteed_to_be_winners.add(var.name())
                        winners = winners.union({var.name()})
                else:
                    small_solutions.append(var)

            len_small_solutions = len(small_solutions)
            len_guides_guaranteed_to_be_winners = len(guides_guaranteed_to_be_winners)
            # --------------------------------------------------------------------------------

            print('{n} guides are guaranteed to be in the final solution (probability of 1.0).'.format(
                n=len_guides_guaranteed_to_be_winners
                )
            )
            print('Exhausting 2^{n} - 1 = {k} possibilities in the worst case.'.format(
                n=len_small_solutions, 
                k=pow(2, len_small_solutions) - 1,
                )
            )

            # ----------------- Combinatorial ---------------------------------
            min_cost = len(self.species)  # Worst case, pick one guide for each species

            for n in range(1, len_small_solutions + 1):
                combos = itertools.combinations(small_solutions, n)

                for combo in combos:
                    possibilities_exhausted += 1

                    I = species_guaranteed_to_be_covered
                    temp_winners = guides_guaranteed_to_be_winners

                    for var in combo:
                        if self.coversets[var.name()][1] - I:
                            I = I.union(self.coversets[var.name()][1])
                            temp_winners = temp_winners.union({var.name()})
                            num_sets = len(combo) + len_guides_guaranteed_to_be_winners

                    if I == self.species and num_sets < min_cost:  # type: ignore
                        min_cost = num_sets  # type: ignore
                        winners = temp_winners
                        break
            
            print('Exhausted {n} possiblities.'.format(n=possibilities_exhausted))
            self.num_exhausted_combos = possibilities_exhausted

            self.exhaustive_time = default_timer() - start_time
            print('\nExhaustive search took {t} seconds'.format(t=round(self.exhaustive_time, 4)))
            # ----------------------------------------------------------------------------

        print('\nSolution set:')
        sol_set: set[int] = set()
        for i, winner in enumerate(winners):
            sol_set = sol_set.union(self.coversets[winner][1])
            score = self.coversets[winner][0]
            self.winner_scores.append(score)

            print('{winner} {score} {c} \t'.format(
                winner=winner, 
                score=score, 
                c=self.coversets[winner][1],
            ).expandtabs(30), end='')

            if i % 2 == 1: print()

        print('\nSize:', len(winners))
        print('Average score: {avg:.2f}'.format(avg=numpy.mean(self.winner_scores)))
        print('Species covered: (union of the solution set): {s}\n'.format(s=sol_set))

        return winners

    
    @property
    def average_score_for_all_trials(self) -> float:
        if self.winner_scores:
            if self.average_score_for_each_trial:
                return numpy.mean(self.average_score_for_each_trial).astype(float)
            else:
                return numpy.mean(self.winner_scores).astype(float)
        else:
            return 0.0
