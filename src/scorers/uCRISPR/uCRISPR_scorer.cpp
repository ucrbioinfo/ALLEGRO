// Find the original paper here:
// uCRISPR - Unified energetics analysis to evaluate the Cas9 on-target
// activity and off-target effects. Version 0.1 Author: Dong Zhang, Travis Hurst,
// Dongsheng Duan & Shi-Jie Chen Email: chenshi@missouri.edu Date: Feb 10, 2019
//
// Find the original open-source project here https://github.com/Vfold-RNA/uCRISPR
//
// Code modified and adapted by Amir to work as a portable library for ALLEGRO.
// No opening and closing files on disk anymore.
// Now directly uses the RNA class from the RNAstructure library.
//
// To build this file, run
// g++ -std=c++17 main.cpp uCRISPR_scorer.cpp RNAstructure/RNA_class/*.o RNAstructure/src/*.o -o uCRISPR_scorer
//
#include <cmath>
#include <mutex>
#include <thread>
#include <stdlib.h>
#include <iostream>
#include <filesystem>

#include "include/uCRISPR_scorer.h"
#include "RNAstructure/RNA_class/RNA.h"

std::mutex mtx; // Mutex to protect shared data.

namespace uCRISPR_scorer
{
    uCRISPR_scorer::uCRISPR_scorer()
    {
        this->ReadParameters();

        // Set RNAstructure data tables environment path.
        std::filesystem::path current_dir = std::filesystem::current_path();
        std::string env_path = current_dir.string() + "/src/scorers/uCRISPR/RNAstructure/data_tables/";
        setenv("DATAPATH", env_path.c_str(), 1);
    }

    uCRISPR_scorer::~uCRISPR_scorer() {}

    // Hard-coded parameters to eliminate the need of disk reads.
    void uCRISPR_scorer::ReadParameters()
    {
        this->parameters_on.insert(std::pair<std::string, double>("AA_-3", 0.0925819));
        this->parameters_on.insert(std::pair<std::string, double>("AC_-3", 0.0922315));
        this->parameters_on.insert(std::pair<std::string, double>("AG_-3", 0.0938087));
        this->parameters_on.insert(std::pair<std::string, double>("AT_-3", 0.0937091));
        this->parameters_on.insert(std::pair<std::string, double>("CA_-3", 0.0949586));
        this->parameters_on.insert(std::pair<std::string, double>("CC_-3", 0.0951256));
        this->parameters_on.insert(std::pair<std::string, double>("CG_-3", 0.0956806));
        this->parameters_on.insert(std::pair<std::string, double>("CT_-3", 0.0921751));
        this->parameters_on.insert(std::pair<std::string, double>("GA_-3", 0.0908114));
        this->parameters_on.insert(std::pair<std::string, double>("GC_-3", 0.09429));
        this->parameters_on.insert(std::pair<std::string, double>("GG_-3", 0.0911885));
        this->parameters_on.insert(std::pair<std::string, double>("GT_-3", 0.0926019));
        this->parameters_on.insert(std::pair<std::string, double>("TA_-3", 0.0903281));
        this->parameters_on.insert(std::pair<std::string, double>("TC_-3", 0.093484));
        this->parameters_on.insert(std::pair<std::string, double>("TG_-3", 0.0922147));
        this->parameters_on.insert(std::pair<std::string, double>("TT_-3", 0.0924238));
        this->parameters_on.insert(std::pair<std::string, double>("AA_-2", 0.0931903));
        this->parameters_on.insert(std::pair<std::string, double>("AC_-2", 0.0920265));
        this->parameters_on.insert(std::pair<std::string, double>("AG_-2", 0.0930193));
        this->parameters_on.insert(std::pair<std::string, double>("AT_-2", 0.0904439));
        this->parameters_on.insert(std::pair<std::string, double>("CA_-2", 0.0918626));
        this->parameters_on.insert(std::pair<std::string, double>("CC_-2", 0.0925978));
        this->parameters_on.insert(std::pair<std::string, double>("CG_-2", 0.0952456));
        this->parameters_on.insert(std::pair<std::string, double>("CT_-2", 0.0954251));
        this->parameters_on.insert(std::pair<std::string, double>("GA_-2", 0.0945837));
        this->parameters_on.insert(std::pair<std::string, double>("GC_-2", 0.0951149));
        this->parameters_on.insert(std::pair<std::string, double>("GG_-2", 0.0931975));
        this->parameters_on.insert(std::pair<std::string, double>("GT_-2", 0.0899964));
        this->parameters_on.insert(std::pair<std::string, double>("TA_-2", 0.0891982));
        this->parameters_on.insert(std::pair<std::string, double>("TC_-2", 0.0938981));
        this->parameters_on.insert(std::pair<std::string, double>("TG_-2", 0.094208));
        this->parameters_on.insert(std::pair<std::string, double>("TT_-2", 0.0936055));
        this->parameters_on.insert(std::pair<std::string, double>("AA_-1", 0.0898541));
        this->parameters_on.insert(std::pair<std::string, double>("AC_-1", 0.093931));
        this->parameters_on.insert(std::pair<std::string, double>("AG_-1", 0.0929465));
        this->parameters_on.insert(std::pair<std::string, double>("AT_-1", 0.0921033));
        this->parameters_on.insert(std::pair<std::string, double>("CA_-1", 0.093045));
        this->parameters_on.insert(std::pair<std::string, double>("CC_-1", 0.0925018));
        this->parameters_on.insert(std::pair<std::string, double>("CG_-1", 0.0944526));
        this->parameters_on.insert(std::pair<std::string, double>("CT_-1", 0.093638));
        this->parameters_on.insert(std::pair<std::string, double>("GA_-1", 0.0941225));
        this->parameters_on.insert(std::pair<std::string, double>("GC_-1", 0.0934688));
        this->parameters_on.insert(std::pair<std::string, double>("GG_-1", 0.0929426));
        this->parameters_on.insert(std::pair<std::string, double>("GT_-1", 0.0951365));
        this->parameters_on.insert(std::pair<std::string, double>("TA_-1", 0.0948378));
        this->parameters_on.insert(std::pair<std::string, double>("TC_-1", 0.0929197));
        this->parameters_on.insert(std::pair<std::string, double>("TG_-1", 0.0909034));
        this->parameters_on.insert(std::pair<std::string, double>("TT_-1", 0.09081));
        this->parameters_on.insert(std::pair<std::string, double>("AA_0", 0.0931161));
        this->parameters_on.insert(std::pair<std::string, double>("AC_0", 0.0907721));
        this->parameters_on.insert(std::pair<std::string, double>("AG_0", 0.0941268));
        this->parameters_on.insert(std::pair<std::string, double>("AT_0", 0.0938444));
        this->parameters_on.insert(std::pair<std::string, double>("CA_0", 0.0923403));
        this->parameters_on.insert(std::pair<std::string, double>("CC_0", 0.091771));
        this->parameters_on.insert(std::pair<std::string, double>("CG_0", 0.0954305));
        this->parameters_on.insert(std::pair<std::string, double>("CT_0", 0.0932796));
        this->parameters_on.insert(std::pair<std::string, double>("GA_0", 0.0942394));
        this->parameters_on.insert(std::pair<std::string, double>("GC_0", 0.0899311));
        this->parameters_on.insert(std::pair<std::string, double>("GG_0", 0.0953627));
        this->parameters_on.insert(std::pair<std::string, double>("GT_0", 0.0917119));
        this->parameters_on.insert(std::pair<std::string, double>("TA_0", 0.0918751));
        this->parameters_on.insert(std::pair<std::string, double>("TC_0", 0.091151));
        this->parameters_on.insert(std::pair<std::string, double>("TG_0", 0.0948152));
        this->parameters_on.insert(std::pair<std::string, double>("TT_0", 0.0938463));
        this->parameters_on.insert(std::pair<std::string, double>("AA_1", 0.0893329));
        this->parameters_on.insert(std::pair<std::string, double>("AC_1", 0.0971859));
        this->parameters_on.insert(std::pair<std::string, double>("AG_1", 0.0931226));
        this->parameters_on.insert(std::pair<std::string, double>("AT_1", 0.0919295));
        this->parameters_on.insert(std::pair<std::string, double>("CA_1", 0.0929137));
        this->parameters_on.insert(std::pair<std::string, double>("CC_1", 0.090097));
        this->parameters_on.insert(std::pair<std::string, double>("CG_1", 0.089183));
        this->parameters_on.insert(std::pair<std::string, double>("CT_1", 0.0914314));
        this->parameters_on.insert(std::pair<std::string, double>("GA_1", 0.0959187));
        this->parameters_on.insert(std::pair<std::string, double>("GC_1", 0.0916874));
        this->parameters_on.insert(std::pair<std::string, double>("GG_1", 0.0942516));
        this->parameters_on.insert(std::pair<std::string, double>("GT_1", 0.0978776));
        this->parameters_on.insert(std::pair<std::string, double>("TA_1", 0.0941265));
        this->parameters_on.insert(std::pair<std::string, double>("TC_1", 0.0962411));
        this->parameters_on.insert(std::pair<std::string, double>("TG_1", 0.0944066));
        this->parameters_on.insert(std::pair<std::string, double>("TT_1", 0.0879081));
        this->parameters_on.insert(std::pair<std::string, double>("AA_2", 0.0893318));
        this->parameters_on.insert(std::pair<std::string, double>("AC_2", 0.0916351));
        this->parameters_on.insert(std::pair<std::string, double>("AG_2", 0.0954458));
        this->parameters_on.insert(std::pair<std::string, double>("AT_2", 0.095879));
        this->parameters_on.insert(std::pair<std::string, double>("CA_2", 0.0951659));
        this->parameters_on.insert(std::pair<std::string, double>("CC_2", 0.092217));
        this->parameters_on.insert(std::pair<std::string, double>("CG_2", 0.0912652));
        this->parameters_on.insert(std::pair<std::string, double>("CT_2", 0.0965632));
        this->parameters_on.insert(std::pair<std::string, double>("GA_2", 0.0945108));
        this->parameters_on.insert(std::pair<std::string, double>("GC_2", 0.0895005));
        this->parameters_on.insert(std::pair<std::string, double>("GG_2", 0.0917777));
        this->parameters_on.insert(std::pair<std::string, double>("GT_2", 0.0951748));
        this->parameters_on.insert(std::pair<std::string, double>("TA_2", 0.0977967));
        this->parameters_on.insert(std::pair<std::string, double>("TC_2", 0.0900949));
        this->parameters_on.insert(std::pair<std::string, double>("TG_2", 0.0934703));
        this->parameters_on.insert(std::pair<std::string, double>("TT_2", 0.0877846));
        this->parameters_on.insert(std::pair<std::string, double>("AA_3", 0.0939707));
        this->parameters_on.insert(std::pair<std::string, double>("AC_3", 0.0934277));
        this->parameters_on.insert(std::pair<std::string, double>("AG_3", 0.0969887));
        this->parameters_on.insert(std::pair<std::string, double>("AT_3", 0.0924182));
        this->parameters_on.insert(std::pair<std::string, double>("CA_3", 0.0929031));
        this->parameters_on.insert(std::pair<std::string, double>("CC_3", 0.0893678));
        this->parameters_on.insert(std::pair<std::string, double>("CG_3", 0.0910223));
        this->parameters_on.insert(std::pair<std::string, double>("CT_3", 0.0901543));
        this->parameters_on.insert(std::pair<std::string, double>("GA_3", 0.0921615));
        this->parameters_on.insert(std::pair<std::string, double>("GC_3", 0.0905665));
        this->parameters_on.insert(std::pair<std::string, double>("GG_3", 0.0941067));
        this->parameters_on.insert(std::pair<std::string, double>("GT_3", 0.0951244));
        this->parameters_on.insert(std::pair<std::string, double>("TA_3", 0.0949984));
        this->parameters_on.insert(std::pair<std::string, double>("TC_3", 0.0929395));
        this->parameters_on.insert(std::pair<std::string, double>("TG_3", 0.0971426));
        this->parameters_on.insert(std::pair<std::string, double>("TT_3", 0.0903212));
        this->parameters_on.insert(std::pair<std::string, double>("AA_4", 0.0916468));
        this->parameters_on.insert(std::pair<std::string, double>("AC_4", 0.0965841));
        this->parameters_on.insert(std::pair<std::string, double>("AG_4", 0.0920573));
        this->parameters_on.insert(std::pair<std::string, double>("AT_4", 0.0937455));
        this->parameters_on.insert(std::pair<std::string, double>("CA_4", 0.0931933));
        this->parameters_on.insert(std::pair<std::string, double>("CC_4", 0.0896192));
        this->parameters_on.insert(std::pair<std::string, double>("CG_4", 0.093097));
        this->parameters_on.insert(std::pair<std::string, double>("CT_4", 0.090392));
        this->parameters_on.insert(std::pair<std::string, double>("GA_4", 0.093783));
        this->parameters_on.insert(std::pair<std::string, double>("GC_4", 0.0937979));
        this->parameters_on.insert(std::pair<std::string, double>("GG_4", 0.0925274));
        this->parameters_on.insert(std::pair<std::string, double>("GT_4", 0.099152));
        this->parameters_on.insert(std::pair<std::string, double>("TA_4", 0.0925122));
        this->parameters_on.insert(std::pair<std::string, double>("TC_4", 0.0957103));
        this->parameters_on.insert(std::pair<std::string, double>("TG_4", 0.0930223));
        this->parameters_on.insert(std::pair<std::string, double>("TT_4", 0.0867732));
        this->parameters_on.insert(std::pair<std::string, double>("AA_5", 0.0897272));
        this->parameters_on.insert(std::pair<std::string, double>("AC_5", 0.0979949));
        this->parameters_on.insert(std::pair<std::string, double>("AG_5", 0.093501));
        this->parameters_on.insert(std::pair<std::string, double>("AT_5", 0.0899123));
        this->parameters_on.insert(std::pair<std::string, double>("CA_5", 0.0938118));
        this->parameters_on.insert(std::pair<std::string, double>("CC_5", 0.0921319));
        this->parameters_on.insert(std::pair<std::string, double>("CG_5", 0.0940882));
        this->parameters_on.insert(std::pair<std::string, double>("CT_5", 0.0956795));
        this->parameters_on.insert(std::pair<std::string, double>("GA_5", 0.0954024));
        this->parameters_on.insert(std::pair<std::string, double>("GC_5", 0.0897733));
        this->parameters_on.insert(std::pair<std::string, double>("GG_5", 0.0942577));
        this->parameters_on.insert(std::pair<std::string, double>("GT_5", 0.0912707));
        this->parameters_on.insert(std::pair<std::string, double>("TA_5", 0.0893433));
        this->parameters_on.insert(std::pair<std::string, double>("TC_5", 0.0950183));
        this->parameters_on.insert(std::pair<std::string, double>("TG_5", 0.0963269));
        this->parameters_on.insert(std::pair<std::string, double>("TT_5", 0.0893742));
        this->parameters_on.insert(std::pair<std::string, double>("AA_6", 0.0877049));
        this->parameters_on.insert(std::pair<std::string, double>("AC_6", 0.0967083));
        this->parameters_on.insert(std::pair<std::string, double>("AG_6", 0.09364));
        this->parameters_on.insert(std::pair<std::string, double>("AT_6", 0.0902315));
        this->parameters_on.insert(std::pair<std::string, double>("CA_6", 0.0951412));
        this->parameters_on.insert(std::pair<std::string, double>("CC_6", 0.0894165));
        this->parameters_on.insert(std::pair<std::string, double>("CG_6", 0.0969462));
        this->parameters_on.insert(std::pair<std::string, double>("CT_6", 0.0934144));
        this->parameters_on.insert(std::pair<std::string, double>("GA_6", 0.0928866));
        this->parameters_on.insert(std::pair<std::string, double>("GC_6", 0.0945298));
        this->parameters_on.insert(std::pair<std::string, double>("GG_6", 0.0946059));
        this->parameters_on.insert(std::pair<std::string, double>("GT_6", 0.0961515));
        this->parameters_on.insert(std::pair<std::string, double>("TA_6", 0.0956894));
        this->parameters_on.insert(std::pair<std::string, double>("TC_6", 0.091851));
        this->parameters_on.insert(std::pair<std::string, double>("TG_6", 0.0936));
        this->parameters_on.insert(std::pair<std::string, double>("TT_6", 0.0850963));
        this->parameters_on.insert(std::pair<std::string, double>("AA_7", 0.0913231));
        this->parameters_on.insert(std::pair<std::string, double>("AC_7", 0.0947777));
        this->parameters_on.insert(std::pair<std::string, double>("AG_7", 0.0931175));
        this->parameters_on.insert(std::pair<std::string, double>("AT_7", 0.0922038));
        this->parameters_on.insert(std::pair<std::string, double>("CA_7", 0.0936672));
        this->parameters_on.insert(std::pair<std::string, double>("CC_7", 0.0901245));
        this->parameters_on.insert(std::pair<std::string, double>("CG_7", 0.0967291));
        this->parameters_on.insert(std::pair<std::string, double>("CT_7", 0.0919848));
        this->parameters_on.insert(std::pair<std::string, double>("GA_7", 0.0938562));
        this->parameters_on.insert(std::pair<std::string, double>("GC_7", 0.0953601));
        this->parameters_on.insert(std::pair<std::string, double>("GG_7", 0.0937199));
        this->parameters_on.insert(std::pair<std::string, double>("GT_7", 0.0958558));
        this->parameters_on.insert(std::pair<std::string, double>("TA_7", 0.0911729));
        this->parameters_on.insert(std::pair<std::string, double>("TC_7", 0.0923409));
        this->parameters_on.insert(std::pair<std::string, double>("TG_7", 0.0915619));
        this->parameters_on.insert(std::pair<std::string, double>("TT_7", 0.089818));
        this->parameters_on.insert(std::pair<std::string, double>("AA_8", 0.090893));
        this->parameters_on.insert(std::pair<std::string, double>("AC_8", 0.0944247));
        this->parameters_on.insert(std::pair<std::string, double>("AG_8", 0.094571));
        this->parameters_on.insert(std::pair<std::string, double>("AT_8", 0.0901307));
        this->parameters_on.insert(std::pair<std::string, double>("CA_8", 0.0930374));
        this->parameters_on.insert(std::pair<std::string, double>("CC_8", 0.0913408));
        this->parameters_on.insert(std::pair<std::string, double>("CG_8", 0.0944616));
        this->parameters_on.insert(std::pair<std::string, double>("CT_8", 0.0937633));
        this->parameters_on.insert(std::pair<std::string, double>("GA_8", 0.0957338));
        this->parameters_on.insert(std::pair<std::string, double>("GC_8", 0.0943693));
        this->parameters_on.insert(std::pair<std::string, double>("GG_8", 0.0885683));
        this->parameters_on.insert(std::pair<std::string, double>("GT_8", 0.0964571));
        this->parameters_on.insert(std::pair<std::string, double>("TA_8", 0.0984968));
        this->parameters_on.insert(std::pair<std::string, double>("TC_8", 0.0902878));
        this->parameters_on.insert(std::pair<std::string, double>("TG_8", 0.0936647));
        this->parameters_on.insert(std::pair<std::string, double>("TT_8", 0.0874132));
        this->parameters_on.insert(std::pair<std::string, double>("AA_9", 0.0936968));
        this->parameters_on.insert(std::pair<std::string, double>("AC_9", 0.0953001));
        this->parameters_on.insert(std::pair<std::string, double>("AG_9", 0.0965613));
        this->parameters_on.insert(std::pair<std::string, double>("AT_9", 0.0926029));
        this->parameters_on.insert(std::pair<std::string, double>("CA_9", 0.0939941));
        this->parameters_on.insert(std::pair<std::string, double>("CC_9", 0.089749));
        this->parameters_on.insert(std::pair<std::string, double>("CG_9", 0.0931962));
        this->parameters_on.insert(std::pair<std::string, double>("CT_9", 0.0934833));
        this->parameters_on.insert(std::pair<std::string, double>("GA_9", 0.0955182));
        this->parameters_on.insert(std::pair<std::string, double>("GC_9", 0.0941275));
        this->parameters_on.insert(std::pair<std::string, double>("GG_9", 0.0902993));
        this->parameters_on.insert(std::pair<std::string, double>("GT_9", 0.0913206));
        this->parameters_on.insert(std::pair<std::string, double>("TA_9", 0.0956822));
        this->parameters_on.insert(std::pair<std::string, double>("TC_9", 0.0930747));
        this->parameters_on.insert(std::pair<std::string, double>("TG_9", 0.090146));
        this->parameters_on.insert(std::pair<std::string, double>("TT_9", 0.0888613));
        this->parameters_on.insert(std::pair<std::string, double>("AA_10", 0.0924818));
        this->parameters_on.insert(std::pair<std::string, double>("AC_10", 0.0933774));
        this->parameters_on.insert(std::pair<std::string, double>("AG_10", 0.0987523));
        this->parameters_on.insert(std::pair<std::string, double>("AT_10", 0.0942799));
        this->parameters_on.insert(std::pair<std::string, double>("CA_10", 0.0928727));
        this->parameters_on.insert(std::pair<std::string, double>("CC_10", 0.0935948));
        this->parameters_on.insert(std::pair<std::string, double>("CG_10", 0.0907724));
        this->parameters_on.insert(std::pair<std::string, double>("CT_10", 0.0950114));
        this->parameters_on.insert(std::pair<std::string, double>("GA_10", 0.0994727));
        this->parameters_on.insert(std::pair<std::string, double>("GC_10", 0.0892946));
        this->parameters_on.insert(std::pair<std::string, double>("GG_10", 0.0906833));
        this->parameters_on.insert(std::pair<std::string, double>("GT_10", 0.0907522));
        this->parameters_on.insert(std::pair<std::string, double>("TA_10", 0.0952541));
        this->parameters_on.insert(std::pair<std::string, double>("TC_10", 0.0925514));
        this->parameters_on.insert(std::pair<std::string, double>("TG_10", 0.0951809));
        this->parameters_on.insert(std::pair<std::string, double>("TT_10", 0.0832817));
        this->parameters_on.insert(std::pair<std::string, double>("AA_11", 0.092561));
        this->parameters_on.insert(std::pair<std::string, double>("AC_11", 0.0953214));
        this->parameters_on.insert(std::pair<std::string, double>("AG_11", 0.09796));
        this->parameters_on.insert(std::pair<std::string, double>("AT_11", 0.0942388));
        this->parameters_on.insert(std::pair<std::string, double>("CA_11", 0.0929185));
        this->parameters_on.insert(std::pair<std::string, double>("CC_11", 0.0896416));
        this->parameters_on.insert(std::pair<std::string, double>("CG_11", 0.0942117));
        this->parameters_on.insert(std::pair<std::string, double>("CT_11", 0.0920464));
        this->parameters_on.insert(std::pair<std::string, double>("GA_11", 0.0990035));
        this->parameters_on.insert(std::pair<std::string, double>("GC_11", 0.0907618));
        this->parameters_on.insert(std::pair<std::string, double>("GG_11", 0.0923331));
        this->parameters_on.insert(std::pair<std::string, double>("GT_11", 0.0932904));
        this->parameters_on.insert(std::pair<std::string, double>("TA_11", 0.0964246));
        this->parameters_on.insert(std::pair<std::string, double>("TC_11", 0.0917348));
        this->parameters_on.insert(std::pair<std::string, double>("TG_11", 0.0897273));
        this->parameters_on.insert(std::pair<std::string, double>("TT_11", 0.0854385));
        this->parameters_on.insert(std::pair<std::string, double>("AA_12", 0.0910947));
        this->parameters_on.insert(std::pair<std::string, double>("AC_12", 0.0969815));
        this->parameters_on.insert(std::pair<std::string, double>("AG_12", 0.0975956));
        this->parameters_on.insert(std::pair<std::string, double>("AT_12", 0.0952359));
        this->parameters_on.insert(std::pair<std::string, double>("CA_12", 0.0910651));
        this->parameters_on.insert(std::pair<std::string, double>("CC_12", 0.0923125));
        this->parameters_on.insert(std::pair<std::string, double>("CG_12", 0.0919526));
        this->parameters_on.insert(std::pair<std::string, double>("CT_12", 0.0921294));
        this->parameters_on.insert(std::pair<std::string, double>("GA_12", 0.0947763));
        this->parameters_on.insert(std::pair<std::string, double>("GC_12", 0.0917749));
        this->parameters_on.insert(std::pair<std::string, double>("GG_12", 0.0913023));
        this->parameters_on.insert(std::pair<std::string, double>("GT_12", 0.0963786));
        this->parameters_on.insert(std::pair<std::string, double>("TA_12", 0.0948079));
        this->parameters_on.insert(std::pair<std::string, double>("TC_12", 0.0936849));
        this->parameters_on.insert(std::pair<std::string, double>("TG_12", 0.0932656));
        this->parameters_on.insert(std::pair<std::string, double>("TT_12", 0.0832557));
        this->parameters_on.insert(std::pair<std::string, double>("AA_13", 0.0928307));
        this->parameters_on.insert(std::pair<std::string, double>("AC_13", 0.0910303));
        this->parameters_on.insert(std::pair<std::string, double>("AG_13", 0.0947911));
        this->parameters_on.insert(std::pair<std::string, double>("AT_13", 0.093092));
        this->parameters_on.insert(std::pair<std::string, double>("CA_13", 0.0924135));
        this->parameters_on.insert(std::pair<std::string, double>("CC_13", 0.092999));
        this->parameters_on.insert(std::pair<std::string, double>("CG_13", 0.088025));
        this->parameters_on.insert(std::pair<std::string, double>("CT_13", 0.101316));
        this->parameters_on.insert(std::pair<std::string, double>("GA_13", 0.0972619));
        this->parameters_on.insert(std::pair<std::string, double>("GC_13", 0.0929587));
        this->parameters_on.insert(std::pair<std::string, double>("GG_13", 0.0892356));
        this->parameters_on.insert(std::pair<std::string, double>("GT_13", 0.09466));
        this->parameters_on.insert(std::pair<std::string, double>("TA_13", 0.097519));
        this->parameters_on.insert(std::pair<std::string, double>("TC_13", 0.0933949));
        this->parameters_on.insert(std::pair<std::string, double>("TG_13", 0.0878261));
        this->parameters_on.insert(std::pair<std::string, double>("TT_13", 0.0882596));
        this->parameters_on.insert(std::pair<std::string, double>("AA_14", 0.0947295));
        this->parameters_on.insert(std::pair<std::string, double>("AC_14", 0.0964921));
        this->parameters_on.insert(std::pair<std::string, double>("AG_14", 0.0947994));
        this->parameters_on.insert(std::pair<std::string, double>("AT_14", 0.0940041));
        this->parameters_on.insert(std::pair<std::string, double>("CA_14", 0.0918242));
        this->parameters_on.insert(std::pair<std::string, double>("CC_14", 0.0903739));
        this->parameters_on.insert(std::pair<std::string, double>("CG_14", 0.0926988));
        this->parameters_on.insert(std::pair<std::string, double>("CT_14", 0.095486));
        this->parameters_on.insert(std::pair<std::string, double>("GA_14", 0.0938239));
        this->parameters_on.insert(std::pair<std::string, double>("GC_14", 0.0881097));
        this->parameters_on.insert(std::pair<std::string, double>("GG_14", 0.0887479));
        this->parameters_on.insert(std::pair<std::string, double>("GT_14", 0.0891962));
        this->parameters_on.insert(std::pair<std::string, double>("TA_14", 0.101033));
        this->parameters_on.insert(std::pair<std::string, double>("TC_14", 0.093477));
        this->parameters_on.insert(std::pair<std::string, double>("TG_14", 0.0928709));
        this->parameters_on.insert(std::pair<std::string, double>("TT_14", 0.0899474));
        this->parameters_on.insert(std::pair<std::string, double>("AA_15", 0.0937616));
        this->parameters_on.insert(std::pair<std::string, double>("AC_15", 0.0999929));
        this->parameters_on.insert(std::pair<std::string, double>("AG_15", 0.0942937));
        this->parameters_on.insert(std::pair<std::string, double>("AT_15", 0.0933619));
        this->parameters_on.insert(std::pair<std::string, double>("CA_15", 0.0945116));
        this->parameters_on.insert(std::pair<std::string, double>("CC_15", 0.0924896));
        this->parameters_on.insert(std::pair<std::string, double>("CG_15", 0.0907899));
        this->parameters_on.insert(std::pair<std::string, double>("CT_15", 0.0906615));
        this->parameters_on.insert(std::pair<std::string, double>("GA_15", 0.0967593));
        this->parameters_on.insert(std::pair<std::string, double>("GC_15", 0.0933444));
        this->parameters_on.insert(std::pair<std::string, double>("GG_15", 0.0876335));
        this->parameters_on.insert(std::pair<std::string, double>("GT_15", 0.0913798));
        this->parameters_on.insert(std::pair<std::string, double>("TA_15", 0.0917906));
        this->parameters_on.insert(std::pair<std::string, double>("TC_15", 0.0975295));
        this->parameters_on.insert(std::pair<std::string, double>("TG_15", 0.0918297));
        this->parameters_on.insert(std::pair<std::string, double>("TT_15", 0.0874839));
        this->parameters_on.insert(std::pair<std::string, double>("AA_16", 0.0918829));
        this->parameters_on.insert(std::pair<std::string, double>("AC_16", 0.0969413));
        this->parameters_on.insert(std::pair<std::string, double>("AG_16", 0.096425));
        this->parameters_on.insert(std::pair<std::string, double>("AT_16", 0.091574));
        this->parameters_on.insert(std::pair<std::string, double>("CA_16", 0.0993696));
        this->parameters_on.insert(std::pair<std::string, double>("CC_16", 0.0983804));
        this->parameters_on.insert(std::pair<std::string, double>("CG_16", 0.092271));
        this->parameters_on.insert(std::pair<std::string, double>("CT_16", 0.0933355));
        this->parameters_on.insert(std::pair<std::string, double>("GA_16", 0.0943085));
        this->parameters_on.insert(std::pair<std::string, double>("GC_16", 0.0926428));
        this->parameters_on.insert(std::pair<std::string, double>("GG_16", 0.085006));
        this->parameters_on.insert(std::pair<std::string, double>("GT_16", 0.0925895));
        this->parameters_on.insert(std::pair<std::string, double>("TA_16", 0.0895233));
        this->parameters_on.insert(std::pair<std::string, double>("TC_16", 0.0911689));
        this->parameters_on.insert(std::pair<std::string, double>("TG_16", 0.0984778));
        this->parameters_on.insert(std::pair<std::string, double>("TT_16", 0.0837171));
        this->parameters_on.insert(std::pair<std::string, double>("AA_17", 0.0872839));
        this->parameters_on.insert(std::pair<std::string, double>("AC_17", 0.104091));
        this->parameters_on.insert(std::pair<std::string, double>("AG_17", 0.0932207));
        this->parameters_on.insert(std::pair<std::string, double>("AT_17", 0.0904883));
        this->parameters_on.insert(std::pair<std::string, double>("CA_17", 0.0927558));
        this->parameters_on.insert(std::pair<std::string, double>("CC_17", 0.100802));
        this->parameters_on.insert(std::pair<std::string, double>("CG_17", 0.0948267));
        this->parameters_on.insert(std::pair<std::string, double>("CT_17", 0.090749));
        this->parameters_on.insert(std::pair<std::string, double>("GA_17", 0.0950984));
        this->parameters_on.insert(std::pair<std::string, double>("GC_17", 0.0923428));
        this->parameters_on.insert(std::pair<std::string, double>("GG_17", 0.0889056));
        this->parameters_on.insert(std::pair<std::string, double>("GT_17", 0.095833));
        this->parameters_on.insert(std::pair<std::string, double>("TA_17", 0.0913612));
        this->parameters_on.insert(std::pair<std::string, double>("TC_17", 0.0932482));
        this->parameters_on.insert(std::pair<std::string, double>("TG_17", 0.101015));
        this->parameters_on.insert(std::pair<std::string, double>("TT_17", 0.0755916));
        this->parameters_on.insert(std::pair<std::string, double>("AA_18", 0.084763));
        this->parameters_on.insert(std::pair<std::string, double>("AC_18", 0.101963));
        this->parameters_on.insert(std::pair<std::string, double>("AG_18", 0.0943017));
        this->parameters_on.insert(std::pair<std::string, double>("AT_18", 0.085472));
        this->parameters_on.insert(std::pair<std::string, double>("CA_18", 0.100182));
        this->parameters_on.insert(std::pair<std::string, double>("CC_18", 0.0972185));
        this->parameters_on.insert(std::pair<std::string, double>("CG_18", 0.0999309));
        this->parameters_on.insert(std::pair<std::string, double>("CT_18", 0.0931527));
        this->parameters_on.insert(std::pair<std::string, double>("GA_18", 0.0960702));
        this->parameters_on.insert(std::pair<std::string, double>("GC_18", 0.0903934));
        this->parameters_on.insert(std::pair<std::string, double>("GG_18", 0.0955619));
        this->parameters_on.insert(std::pair<std::string, double>("GT_18", 0.0959427));
        this->parameters_on.insert(std::pair<std::string, double>("TA_18", 0.0869178));
        this->parameters_on.insert(std::pair<std::string, double>("TC_18", 0.0900774));
        this->parameters_on.insert(std::pair<std::string, double>("TG_18", 0.0992909));
        this->parameters_on.insert(std::pair<std::string, double>("TT_18", 0.0763758));
        this->parameters_on.insert(std::pair<std::string, double>("AA_19", 0.0867305));
        this->parameters_on.insert(std::pair<std::string, double>("AC_19", 0.0919618));
        this->parameters_on.insert(std::pair<std::string, double>("AG_19", 0.100213));
        this->parameters_on.insert(std::pair<std::string, double>("AT_19", 0.089028));
        this->parameters_on.insert(std::pair<std::string, double>("CA_19", 0.0991234));
        this->parameters_on.insert(std::pair<std::string, double>("CC_19", 0.0914674));
        this->parameters_on.insert(std::pair<std::string, double>("CG_19", 0.102476));
        this->parameters_on.insert(std::pair<std::string, double>("CT_19", 0.0865849));
        this->parameters_on.insert(std::pair<std::string, double>("GA_19", 0.0994937));
        this->parameters_on.insert(std::pair<std::string, double>("GC_19", 0.090148));
        this->parameters_on.insert(std::pair<std::string, double>("GG_19", 0.0978674));
        this->parameters_on.insert(std::pair<std::string, double>("GT_19", 0.101576));
        this->parameters_on.insert(std::pair<std::string, double>("TA_19", 0.0887301));
        this->parameters_on.insert(std::pair<std::string, double>("TC_19", 0.0828263));
        this->parameters_on.insert(std::pair<std::string, double>("TG_19", 0.102761));
        this->parameters_on.insert(std::pair<std::string, double>("TT_19", 0.0766263));
        this->parameters_on.insert(std::pair<std::string, double>("AA_20", 0.0908818));
        this->parameters_on.insert(std::pair<std::string, double>("AC_20", 0.0930678));
        this->parameters_on.insert(std::pair<std::string, double>("AG_20", 0.0979645));
        this->parameters_on.insert(std::pair<std::string, double>("AT_20", 0.0921636));
        this->parameters_on.insert(std::pair<std::string, double>("CA_20", 0.0872556));
        this->parameters_on.insert(std::pair<std::string, double>("CC_20", 0.092153));
        this->parameters_on.insert(std::pair<std::string, double>("CG_20", 0.0925016));
        this->parameters_on.insert(std::pair<std::string, double>("CT_20", 0.0844932));
        this->parameters_on.insert(std::pair<std::string, double>("GA_20", 0.0996894));
        this->parameters_on.insert(std::pair<std::string, double>("GC_20", 0.100152));
        this->parameters_on.insert(std::pair<std::string, double>("GG_20", 0.10182));
        this->parameters_on.insert(std::pair<std::string, double>("GT_20", 0.101655));
        this->parameters_on.insert(std::pair<std::string, double>("TA_20", 0.0940063));
        this->parameters_on.insert(std::pair<std::string, double>("TC_20", 0.0862813));
        this->parameters_on.insert(std::pair<std::string, double>("TG_20", 0.0909051));
        this->parameters_on.insert(std::pair<std::string, double>("TT_20", 0.082623));
        this->parameters_on.insert(std::pair<std::string, double>("AA_21", 0.0938151));
        this->parameters_on.insert(std::pair<std::string, double>("AC_21", 0.0930845));
        this->parameters_on.insert(std::pair<std::string, double>("AG_21", 0.0914588));
        this->parameters_on.insert(std::pair<std::string, double>("AT_21", 0.0934747));
        this->parameters_on.insert(std::pair<std::string, double>("CA_21", 0.0921513));
        this->parameters_on.insert(std::pair<std::string, double>("CC_21", 0.0954392));
        this->parameters_on.insert(std::pair<std::string, double>("CG_21", 0.0894332));
        this->parameters_on.insert(std::pair<std::string, double>("CT_21", 0.0946306));
        this->parameters_on.insert(std::pair<std::string, double>("GA_21", 0.097474));
        this->parameters_on.insert(std::pair<std::string, double>("GC_21", 0.0957932));
        this->parameters_on.insert(std::pair<std::string, double>("GG_21", 0.0912505));
        this->parameters_on.insert(std::pair<std::string, double>("GT_21", 0.0986732));
        this->parameters_on.insert(std::pair<std::string, double>("TA_21", 0.0880918));
        this->parameters_on.insert(std::pair<std::string, double>("TC_21", 0.0939867));
        this->parameters_on.insert(std::pair<std::string, double>("TG_21", 0.0840529));
        this->parameters_on.insert(std::pair<std::string, double>("TT_21", 0.0948037));
        this->parameters_on.insert(std::pair<std::string, double>("AA_22", 0.0888127));
        this->parameters_on.insert(std::pair<std::string, double>("AC_22", 0.0954234));
        this->parameters_on.insert(std::pair<std::string, double>("AG_22", 0.0943625));
        this->parameters_on.insert(std::pair<std::string, double>("AT_22", 0.0929336));
        this->parameters_on.insert(std::pair<std::string, double>("CA_22", 0.0973499));
        this->parameters_on.insert(std::pair<std::string, double>("CC_22", 0.0903041));
        this->parameters_on.insert(std::pair<std::string, double>("CG_22", 0.0956434));
        this->parameters_on.insert(std::pair<std::string, double>("CT_22", 0.0950063));
        this->parameters_on.insert(std::pair<std::string, double>("GA_22", 0.0893314));
        this->parameters_on.insert(std::pair<std::string, double>("GC_22", 0.0897983));
        this->parameters_on.insert(std::pair<std::string, double>("GG_22", 0.0888341));
        this->parameters_on.insert(std::pair<std::string, double>("GT_22", 0.0882316));
        this->parameters_on.insert(std::pair<std::string, double>("TA_22", 0.09617));
        this->parameters_on.insert(std::pair<std::string, double>("TC_22", 0.097431));
        this->parameters_on.insert(std::pair<std::string, double>("TG_22", 0.0949289));
        this->parameters_on.insert(std::pair<std::string, double>("TT_22", 0.0930524));
        this->parameters_on.insert(std::pair<std::string, double>("AA_23", 0.0900106));
        this->parameters_on.insert(std::pair<std::string, double>("AC_23", 0.0943862));
        this->parameters_on.insert(std::pair<std::string, double>("AG_23", 0.0967805));
        this->parameters_on.insert(std::pair<std::string, double>("AT_23", 0.0904866));
        this->parameters_on.insert(std::pair<std::string, double>("CA_23", 0.0949282));
        this->parameters_on.insert(std::pair<std::string, double>("CC_23", 0.0903805));
        this->parameters_on.insert(std::pair<std::string, double>("CG_23", 0.0972862));
        this->parameters_on.insert(std::pair<std::string, double>("CT_23", 0.090362));
        this->parameters_on.insert(std::pair<std::string, double>("GA_23", 0.0941261));
        this->parameters_on.insert(std::pair<std::string, double>("GC_23", 0.0940775));
        this->parameters_on.insert(std::pair<std::string, double>("GG_23", 0.092332));
        this->parameters_on.insert(std::pair<std::string, double>("GT_23", 0.0932333));
        this->parameters_on.insert(std::pair<std::string, double>("TA_23", 0.0939149));
        this->parameters_on.insert(std::pair<std::string, double>("TC_23", 0.0910084));
        this->parameters_on.insert(std::pair<std::string, double>("TG_23", 0.0968531));
        this->parameters_on.insert(std::pair<std::string, double>("TT_23", 0.0874474));
        this->parameters_on.insert(std::pair<std::string, double>("wdG", 0.00116983));
    }

    double uCRISPR_scorer::GetSgRNAEnergy(const std::string &sequence)
    {
        std::string sline = "NNNNNNNNNNNNNNNNNNNNGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCU";
        std::string sgRNA = sline.substr(20);
        std::string gRNA = sequence;

        while (gRNA.find("T") != std::string::npos)
        {
            int found = gRNA.find("T");
            gRNA.replace(found, 1, "U");
        }

        sgRNA = gRNA + sgRNA;

        RNA *strand = new RNA(sgRNA.c_str());
        strand->PartitionFunction();
        double energy = strand->GetEnsembleEnergy();

        return energy;
    }

    double uCRISPR_scorer::score_guide(const std::string &sequence)
    {
        std::string dna = "";
        std::string rna = "";
        std::string pam = "";

        if (sequence.size() == 23)
        {
            dna = sequence;
            rna = sequence.substr(0, 20);
            pam = sequence.substr(21);
        }
        else if (sequence.size() == 30)
        {
            dna = sequence.substr(0, 25);
            dna += sequence.substr(27);
            rna = sequence.substr(4, 20);
            pam = sequence.substr(25, 2);
        }

        double E_C = 0.0;
        if (sequence.size() == 23)
        {
            for (int ii = 0; ii < 20; ii++)
            {
                std::string pname = dna.substr(ii, 2) + "_" + std::to_string(ii + 1);
                E_C += parameters_on[pname];
            }
        }
        else if (sequence.size() == 30)
        {
            for (int ii = -3; ii < 24; ii++)
            {
                std::string pname = dna.substr(ii + 3, 2) + "_" + std::to_string(ii);
                E_C += parameters_on[pname];
            }
        }

        double EsgRNA = GetSgRNAEnergy(rna);
        EsgRNA *= parameters_on["wdG"];

        double Epam = 1.0;
        if (pam == "GG")
        {
            Epam = 1.0;
        }
        else if (pam == "AG")
        {
            Epam = 0.259259259;
        }
        else if (pam == "CG")
        {
            Epam = 0.107142857;
        }
        else if (pam == "GA")
        {
            Epam = 0.069444444;
        }
        else if (pam == "GC")
        {
            Epam = 0.022222222;
        }
        else if (pam == "GT")
        {
            Epam = 0.016129032;
        }
        else if (pam == "TG")
        {
            Epam = 0.038961039;
        }
        else
        {
            Epam = 0.0;
        }

        return std::exp(E_C + EsgRNA) * Epam;
    }

    // Thread function to calculate scores for a range of sequences.
    void uCRISPR_scorer::threaded_scorer(const std::vector<std::string> &sequences, std::vector<double> &scores, size_t start, size_t end)
    {
        for (size_t i = start; i < end; ++i)
        {
            double score = this->score_guide(sequences[i]);

            std::lock_guard<std::mutex> lock(mtx); // Lock the mutex to protect shared data.
            scores[i] = score;                     // Store the score in the corresponding index.
        }
    }

    void uCRISPR_scorer::score_guides(const std::vector<std::string> &sequences)
    {
        std::vector<double> scores(sequences.size());

        // Get the available number of threads to use.
        size_t num_threads = std::thread::hardware_concurrency();

        std::vector<std::thread> threads;
        threads.reserve(num_threads);

        // Calculate chunk size for each thread.
        size_t chunk_size = sequences.size() / num_threads;

        // Distribute the work among threads.
        for (size_t i = 0; i < num_threads; ++i)
        {
            size_t start = i * chunk_size;
            size_t end = (i == num_threads - 1) ? sequences.size() : start + chunk_size;

            threads.emplace_back(&uCRISPR_scorer::threaded_scorer, this, std::ref(sequences), std::ref(scores), start, end);
        }

        for (auto &thread : threads)
        {
            thread.join();
        }

        // Print scores.
        for (size_t i = 0; i < sequences.size(); ++i)
        {
            std::cout << scores[i] << std::endl;
        }
    }
}