This is a tailored version of uCRISPR. More details are in uCRISPR_scorer.cpp

Any changes to any of these .cpp and/or .pyx files must be compiled before ALLEGRO can see them.
After you make your changes, run:

```
g++ -std=c++17 main.cpp uCRISPR_scorer.cpp RNAstructure/RNA_class/*.o RNAstructure/src/*.o -o uCRISPR_scorer