#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "./comet_input.h"
#include "./cx_input.h"


using namespace std;

//Energy range for chi squared calculations
double chi_energy_start = 0.370;
double chi_energy_end = 0.95;

int chi_squared(double energy_start, double energy_end, double energy_step, int comet_number){

    //inputs comet name from comet_input.h
    string comet_name = input_comet_name[comet_number - 1];

    //calculates numbers of steps and creates appropriately-sized array
    double energy = energy_start;
    int energy_row = (energy_end - energy_start)/energy_step + 1;

    //input observation file
    string observation_name = "../Inputs/Observations/" + comet_name + "_intensity_Chandra_yerr.dat";
    ifstream observation_file(observation_name.c_str());

    int obs_row;
    observation_file >> obs_row;

    double input_obs[obs_row][3];
    double obs_energy[obs_row];
    double obs_intensity[obs_row];
    double obs_intensity_err[obs_row];
    for( int i=0; i<obs_row; i++ ){
        for( int j=0; j<3; j++ ){
            observation_file >> input_obs[i][j]; }
        obs_energy[i] = input_obs[i][0];
        obs_intensity[i] = input_obs[i][1];
        obs_intensity_err[i] = input_obs[i][2];

    }
	observation_file.close();

    //input model file
    string input_name_model = "../Results/" + comet_name + "/ACIS_CX_spectrum_"+ comet_name + ".dat";
    ifstream model_file(input_name_model.c_str());

    double input_model[energy_row][2];
    double model_energy[energy_row];
    double model_intensity[energy_row];
    for( int i=0; i<energy_row; i++ ){
        for( int j=0; j<2; j++ ){
            model_file >> input_model[i][j]; }
        model_energy[i] = input_model[i][0];
        model_intensity[i] = input_model[i][1];
    }
    model_file.close();

    gsl_interp_accel *accel_ptr_model;
    gsl_spline *spline_ptr_model;
    accel_ptr_model = gsl_interp_accel_alloc ();
    spline_ptr_model = gsl_spline_alloc (gsl_interp_cspline, energy_row);
    gsl_spline_init (spline_ptr_model, model_energy, model_intensity, energy_row);


    //Calculate total number of rows that fall within chi-squared energy range
    int chi_row = 0;
    for ( int i=0; i<obs_row; i++ ){
        if (chi_energy_start <= obs_energy[i] && obs_energy[i] <= chi_energy_end) chi_row++; }

    //Define chi squared variables
    double chi_energy[chi_row];
    double chi_squared[chi_row];
    double delta_chi_squared[chi_row];;
    double chi_squared_total = 0;

    //Calculate chi squared and chi squared residual

    int row_temp = 0;
    for ( int i = 0; i<obs_row; i++ ){
        if ( chi_energy_start <= obs_energy[i] && obs_energy[i] <= chi_energy_end ){
            chi_energy[row_temp] = obs_energy[i];
            double model_temp = gsl_spline_eval (spline_ptr_model, obs_energy[i], accel_ptr_model);

            //chi_squared[row_temp] = (obs_intensity[i] - model_temp)*(obs_intensity[i] - model_temp)/model_temp;
            //chi_squared_total += (obs_intensity[i] - model_temp)*(obs_intensity[i] - model_temp)/model_temp;
            chi_squared[row_temp] = (obs_intensity[i] - model_temp)*(obs_intensity[i] - model_temp)/(obs_intensity_err[i] * obs_intensity_err[i]);
            chi_squared_total += chi_squared[row_temp];

            if ( row_temp > 0 ) delta_chi_squared[row_temp] = chi_squared[row_temp] - chi_squared[row_temp-1];
            row_temp++;
        }
    }

    cout << comet_name << " chi squared: " << chi_squared_total/(row_temp - cx_variables) << endl;

    //Outputs the chi squared and chi squared residual arrays
    string output_name_chi_res = "../Results/" + comet_name + "/Chandra_CX_chi_squared_residual_"+ comet_name + ".dat";
    ofstream output_file_chi_res(output_name_chi_res.c_str());
    string output_name_chi = "../Results/" + comet_name + "/Chandra_CX_chi_squared_"+ comet_name + ".dat";
    ofstream output_file_chi(output_name_chi.c_str());
    for( int i=0; i<chi_row; i++ ){
        output_file_chi_res << scientific << chi_energy[i] << " " << delta_chi_squared[i] << endl;
        output_file_chi << scientific << chi_energy[i] << " " << chi_squared[i] << endl;
    }
    output_file_chi_res.close();
    output_file_chi.close();


    gsl_spline_free (spline_ptr_model);
    gsl_interp_accel_free (accel_ptr_model);

    return 0;
}
