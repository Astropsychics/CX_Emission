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
    int model_number = 0;
    int sw_speed = 0;
    bool chi_square_bool = 0;

    int ISON_number = 0;
    int Tempel1_number = 0;

	//Prompt user to select which model to use
    cout << "\nWelcome to the Cometary CX Model.\n\n"
         << "Please select which CX model you would like to use:\n"
         << "0.Kharchenko Group Model      1.Stancil Group Model\n";
    cin >> model_number;
    cout << endl;
    if ( 0 > model_number || model_number > 1 ){
		cout << "I'm sorry, Dave. I can't do that.\n";
		return 0; }

    //If Stancil Model is selected, user inputs desired SW speed
    if ( model_number == 1 ){
        cout << "Please select your preferred SW speed:\n"
             << "1.100km/s  2.200km/s  3.300km/s  4.400km/s  5.500km/s\n"
             << "6.600km/s  7.700km/s  8.800km/s  9.900km/s  10.1000km/s\n";
        cin >> sw_speed;
        sw_speed = sw_speed - 1;
        cout << endl;
    }


    //Prompt user to select which comet to analyze
    cout << "Which comet would you like to analyze?\n"
         << "1.8P  2.Encke  3.IZ  4.LS4  5.MH\n"
         << "6.PanSTARRS  7.ISON  8.Tempel 1\n";
	cin >> comet_number;
    cout << endl;
	if ( 1 > comet_number || comet_number > 8 ){
		cout << "I'm sorry, Dave. I can't do that.\n";
		return 0; }

    //Prompts user to select specific Tempel 1 observation
    if ( comet_number == 8 ){
        cout << "Select which Tempel 1 Visit you want to model:\n"
             << "1.Average  2.Visit 1  3.Visit 2  4.Visit 3\n"
             << "5.Visit 4  6.Visit 5  7.Visit 6  8.Visit 7\n";
        cin >> Tempel1_number;
        cout << endl;
        comet_number += Tempel1_number + 2;
    }

    //Prompts user to select specific ISON observation
    if ( comet_number == 7 ){
        cout << "Select which ISON Visit you want to model:\n"
             << "1.Average  2.Visit 1  3.Visit 2  4.Visit 3\n";
        cin >> ISON_number;
        cout << endl;
        comet_number += ISON_number - 1;
    }

    //Prompts user to calculate chi-square for all comets with the necessary data
    if ( comet_number >= 6 && comet_number <= 18 ){
        cout << "Good news, everyone!\n"
             << "Chi-square can be calculated for the comet you selected.\n"
             << "For these calculations to be performed, please enter '1'.\n"
             << "Otherwise, enter '0'.\n";
        cin >> chi_square_bool;
        cout << endl;
    }

    //performs cx calculations for user-selected model. The resulting model
    //will be scaled to the comet later
    cx_calculations(energy_start, energy_end, energy_step, width, comet_number,
                    model_number, sw_speed);

    //calculates scaling factor for the comet
    cx_compare(energy_start, energy_end, energy_step, comet_number);

    //calculates modeled ACIS and HRC spectra for the comet
    acis_calculations(energy_start, energy_end, energy_step, comet_number);
    hrc_calculations(energy_start, energy_end, energy_step, comet_number);

    //If selected by the user, calculates chi-square
    if ( chi_square_bool == 1 ){
        chi_squared(energy_start, energy_end, energy_step, comet_number, model_number); }

    return 0;
}
