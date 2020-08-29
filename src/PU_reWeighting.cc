#include "PU_reWeighting.h"

PU::PU(){}
PU::~PU(){}

double PU::PU_2016(int NumTrueInteraction, TString type_str){
if (NumTrueInteraction < 0 || NumTrueInteraction > 74 ) {return 1;}
else {
	if (type_str == "nominal") {return pu2016_nominal[NumTrueInteraction];}
	else if (type_str == "up") {return pu2016_up[NumTrueInteraction];}
	else if (type_str == "down") {return pu2016_down[NumTrueInteraction];}
	else {std::cout<<"Error pu string!"<<std::endl; return 1.0;}
}
}
		
double PU::PU_2017(int NumTrueInteraction, TString type_str){
if (NumTrueInteraction < 0 || NumTrueInteraction > 98 ) {return 1;}
else {
        if (type_str == "nominal") {return pu2017_nominal[NumTrueInteraction];}
        else if (type_str == "up") {return pu2017_up[NumTrueInteraction];}
        else if (type_str == "down") {return pu2017_down[NumTrueInteraction];}
	else {std::cout<<"Error pu string!"<<std::endl; return 1.0;}
}
}
		
double PU::PU_2018(int NumTrueInteraction, TString type_str){
if (NumTrueInteraction < 0 || NumTrueInteraction > 98 ) {return 1;}
else {
        if (type_str == "nominal") {return pu2018_nominal[NumTrueInteraction];}
        else if (type_str == "up") {return pu2018_up[NumTrueInteraction];}
        else if (type_str == "down") {return pu2018_down[NumTrueInteraction];}
	else {std::cout<<"Error pu string!"<<std::endl; return 1.0;}
}
}
				
