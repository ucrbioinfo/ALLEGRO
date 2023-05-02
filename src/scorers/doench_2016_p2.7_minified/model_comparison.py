import numpy
import pandas
import featurization


def concatenate_feature_sets(feature_sets):
    '''
    Given a dictionary of sets of features, each in a Pandas.DataFrame,
    concatenate them together to form one big np.array, and get the dimension
    of each set
    Returns: inputs, dim
    '''

    N = feature_sets[feature_sets.keys()[0]].shape[0]
    inputs = numpy.zeros((N, 0))
    for set in feature_sets.keys():
        inputs_set = feature_sets[set].values
        inputs = numpy.hstack((inputs, inputs_set))

    return inputs


def predict(seqs, model):
    model, learn_options = model

    learn_options["V"] = 2

    Xdf = pandas.DataFrame({
        '30mer': seqs,
        'Strand': ['NA'] * len(seqs)
    })

    feature_sets = featurization.featurize_data(Xdf, learn_options)
    inputs = concatenate_feature_sets(feature_sets)

    # call to scikit-learn, returns a vector of predicted values
    return model.predict(inputs)