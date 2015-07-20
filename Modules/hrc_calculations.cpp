#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <vector>

#include "./comet_input.h"

using namespace std;

int hrc_calculations(double energy_start, double energy_end, double energy_step, int comet_number){

    //inputs comet name from comet_input.h
    string comet_name = input_comet_name[comet_number - 1];

   //calculates numbers of steps and creates appropriately-sized array
   double energy = energy_start;
   int energy_row = (energy_end - energy_start)/energy_step + 1;


   vector<double> input_energy(energy_row, 0);

   for ( int x = 0; x < energy_row; x++ ){
       input_energy[x] = energy;
       energy += energy_step; }

   //input effective area
   fstream input_area_file("../Inputs/Effective_Area/hrc_effective_area.dat", fstream::in);

   int area_row;
   input_area_file >> area_row;

   double input_area[area_row][2];
   for( int i=0; i<area_row; i++ ){
       for( int j=0; j<2; j++ ){
           input_area_file >> input_area[i][j]; } }
   input_area_file.close();


   //input comet CX spectrum
   string input_spectrum_file_name = "../Results/" + comet_name + "/CX_spectrum_"+ comet_name + ".dat";
   ifstream input_spectrum_file(input_spectrum_file_name.c_str());

   double input_CX[energy_row][2];
   for( int i=0; i<energy_row; i++ ){
       for( int j=0; j<2; j++ ){
           input_spectrum_file >> input_CX[i][j]; } }
   input_spectrum_file.close();


   //outputs re-scaled HRC response spectrum
   string output_name_Chandra = "../Results/" + comet_name + "/HRC_CX_response_spectrum_"+ comet_name + ".dat";
   ofstream output_file_Chandra(output_name_Chandra.c_str());
   for( int i=0; i<energy_row; i++ ){
       int z = 0;
       while ( input_area[z][1] <= input_energy[i] ) z++;

       output_file_Chandra << scientific << input_CX[i][0] << " " << input_CX[i][1] * input_area[z][1] << endl;
   }
   output_file_Chandra.close();


   return 0;
}
