#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "./comet_input.h"

using namespace std;

int cx_compare(double energy_start, double energy_end, double energy_step, int comet_number){

    //inputs comet name from comet_input.h
    string comet_name = input_comet_name[comet_number - 1];

    //calculates numbers of steps and creates appropriately-sized array
    double energy = energy_start;
    int energy_row = (energy_end - energy_start)/energy_step + 1;

    //input observation file
    string observation_name = "../Inputs/Observations/" + comet_name + "_intensity_Chandra.dat";
    ifstream observation_file(observation_name.c_str());

    int obs_row;
    observation_file >> obs_row;

    double input_obs[obs_row][2];
    double obs_energy[obs_row];
    double obs_intensity[obs_row];
    for( int i=0; i<obs_row; i++ ){
        for( int j=0; j<2; j++ ){
            observation_file >> input_obs[i][j]; }
        obs_energy[i] = input_obs[i][0];
        obs_intensity[i] = input_obs[i][1];
    }
	observation_file.close();

    gsl_interp_accel *accel_ptr_obs;
    gsl_spline *spline_ptr_obs;
    accel_ptr_obs = gsl_interp_accel_alloc ();
    spline_ptr_obs = gsl_spline_alloc (gsl_interp_cspline, obs_row);
    gsl_spline_init (spline_ptr_obs, obs_energy, obs_intensity, obs_row);


    //calculates numbers of steps and creates appropriately-sized array
    vector<double> input_energy(energy_row, 0);

    for ( int x = 0; x < energy_row; x++ ){
        input_energy[x] = energy;
        energy += energy_step; }

    //input effective area
    fstream input_area_file("../Inputs/Effective_Area/acis_effective_area.dat", fstream::in);

    int area_row;
    input_area_file >> area_row;

    double input_area[area_row][3];
    for( int i=0; i<area_row; i++ ){
        for( int j=0; j<3; j++ ){
            input_area_file >> input_area[i][j]; } }
    input_area_file.close();

    //inputs CX model
    fstream model_file("../Results/CX_spectrum.dat", fstream::in);

    double input_CX[energy_row][2];
    double spectrum_energy[energy_row];
    double spectrum_intensity[energy_row];

    for( int i=0; i<energy_row; i++ ){
        for( int j=0; j<2; j++ ){
            model_file >> input_CX[i][j]; }
        spectrum_energy[i] = input_CX[i][0];
        spectrum_intensity[i] = input_CX[i][1];}
    model_file.close();


    //multiplies model with ACIS effective area
    double model_intensity[energy_row];

    for (int i=0; i<energy_row; i++){
        int j = 0;
        while ( input_area[j][1] <= input_energy[i] ) j++;

        model_intensity[i] = spectrum_intensity[i] * input_area[j][2];
    }

    gsl_interp_accel *accel_ptr_model;
    gsl_spline *spline_ptr_model;
    accel_ptr_model = gsl_interp_accel_alloc ();
    spline_ptr_model = gsl_spline_alloc (gsl_interp_cspline, energy_row);
    gsl_spline_init (spline_ptr_model, spectrum_energy, model_intensity, energy_row);


    //Calculates the average ratio between the observational data and the model to find an
    //optimal scaling factor. Scaling is broken into 'low' and 'high'ranges to
    //account for higher confidence in the model at lower energies due to a larger
    //abundance of modeled emission lines
    double scaling_factor_low = 0;
    double counter_low = 0;
    double scaling_factor_high = 0;
    double counter_high = 0;

    for ( float energy = 0.400 ; energy <= 0.700; energy += 0.02 ){
        double obs_temp = gsl_spline_eval (spline_ptr_obs, energy, accel_ptr_obs);
        double model_temp = gsl_spline_eval (spline_ptr_model, energy, accel_ptr_model);

        scaling_factor_low += obs_temp / model_temp;
        counter_low++;
    }
    for ( float energy = 0.700 ; energy <= 0.900; energy += 0.02 ){
        double obs_temp = gsl_spline_eval (spline_ptr_obs, energy, accel_ptr_obs);
        double model_temp = gsl_spline_eval (spline_ptr_model, energy, accel_ptr_model);

        scaling_factor_high += obs_temp / model_temp;
        counter_high++;
    }

    //Performs weighted average on the two scaling factors to account for higher
    //confidence in lower energies; currently set at 80/20 ratio due to ratio
    //of emission lines present within each energy region
    double average_scaling_factor = 0.8*scaling_factor_low/counter_low +
        0.2*scaling_factor_high/counter_high;

    //outputs re-scaled intensity spectrum
    string output_name_spectrum = "../Results/" + comet_name + "/CX_spectrum_"+ comet_name + ".dat";
    ofstream output_file_spectrum(output_name_spectrum.c_str());
    for( int i=0; i<energy_row; i++ ){
        output_file_spectrum << scientific << spectrum_energy[i] << " "
            << spectrum_intensity[i]*average_scaling_factor << endl;
    }
    output_file_spectrum.close();

    gsl_spline_free (spline_ptr_obs);
    gsl_interp_accel_free (accel_ptr_obs);

    gsl_spline_free (spline_ptr_model);
    gsl_interp_accel_free (accel_ptr_model);

    return 0;
}
