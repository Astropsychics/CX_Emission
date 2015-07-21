#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>

#include "./acis_calculations.cpp"
#include "./chi_squared.cpp"
#include "./compare.cpp"
#include "./cx_calculations.cpp"
#include "./hrc_calculations.cpp"


using namespace std;

//defines energy range and step size
double energy_start = 0.10;
double energy_end = 2.0;
double energy_step = 0.01;
double width = 0.050;

int main(int argc, char *argv[]){

    int comet_number = 0;
    bool chi_square_bool = 0;

	//Prompt user to select which comet to analyze
    cout << "\n Welcome to the Cometary CX Model. \n \n"
         << "Please select the number of which comet you'd like to analyze. \n"
         << "1.8P  2.Encke  3.IZ  4.LS4 \n"
         << "5.MH  6.ISON  7.PanSTARRS \n";
	cin >> comet_number;
    cout << endl;

	if ( 1 > comet_number || comet_number > 7 ){
		cout << "I'm sorry, Dave. I can't do that. \n";
		return 0; }

    //Prompts user to calculate chi-square for all comets with the necessary data
    if ( comet_number == 6 || comet_number == 7 ){
        cout << "Good news, everyone! \n"
             << "Chi-square can be calculated for the comet you selected. \n"
             << "For these calculations to be performed, please enter '1'. \n"
             << "Otherwise, enter '0'. \n";
        cin >> chi_square_bool;
        cout << endl;
    }

    //performs cx calculations that will be scaled to the comet later
    cx_calculations(energy_start, energy_end, energy_step, width, comet_number);

    //calculates scaling factor for the comet
    cx_compare(energy_start, energy_end, energy_step, comet_number);

    //calculates modeled ACIS and HRC spectra for the comet
    acis_calculations(energy_start, energy_end, energy_step, comet_number);
    hrc_calculations(energy_start, energy_end, energy_step, comet_number);

    //If selected by the user, calculates chi-square
    if ( chi_square_bool == 1 ){
        chi_squared(energy_start, energy_end, energy_step, comet_number); }

    return 0;
}
