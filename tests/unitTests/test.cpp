/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "test.h"

using namespace mat_fun;
using namespace std;


TEST(UnitTestIso_1, nHK) {
    // Step 1: define parameters
    auto matType = consts::ConstitutiveModelType::stIso_nHook;   // Material_model: options refer to consts.h 
    auto volType = consts::ConstitutiveModelType::stVol_ST91;   // Dilational_penalty_model
    double E = 240.56596e6;   // Elasticity_modulus
    double nu = 0.5;   // Poisson_ratio
    double pen = 4e9;   // Penalty_parameter
    double C01;   // additional parameter to C10 (optional)

    // Step 2: construct test object
    UnitTestIso nHK(matType, E, nu, volType, pen, C01);

    // Step 3: define the input 
    double F[3][3] = {};
    F[0][0] = 1.0; F[1][1] = 1.0; F[2][2] = 1.0;   // set to Identity

    // Step 4: define the reference output 
    double S_ref[3][3] = {};
    double Dm_ref[6][6] = {};

    // Step 5: run unit test
    nHK.runUnitTest(F, S_ref, Dm_ref);
      
}

TEST(UnitTestIso_2, MR) {
    // Step 1: define parameters
    auto matType = consts::ConstitutiveModelType::stIso_MR;   // Material_model: options refer to consts.h 
    auto volType = consts::ConstitutiveModelType::stVol_ST91;   // Dilational_penalty_model
    double E = 1e6;   // Elasticity_modulus
    double nu = 0.495;   // Poisson_ratio
    double pen = 4e9;   // Penalty_parameter
    double C01 = 0.1;   // additional parameter to C10 (optional)

    // Step 2: construct test object
    UnitTestIso MR(matType, E, nu, volType, pen, C01);

    // Step 3: define the input 
    double F[3][3] = {};
    F[0][0] = 1.0; F[1][1] = 1.0; F[2][2] = 1.0;   // set to Identity

    // Step 4: define the reference output 
    double S_ref[3][3] = {};
    double Dm_ref[6][6] = {};

    // Step 5: run unit test
    MR.runUnitTest(F, S_ref, Dm_ref);
      
}

TEST(UnitTestIso_3, CANN) {
    // Step 1: define parameters
    std::cout << 'test'<<std::ends;
    auto matType = consts::ConstitutiveModelType::stAnisoHyper_Inv;   // Material_model: options refer to consts.h 
    auto volType = consts::ConstitutiveModelType::stVol_ST91;   // Dilational_penalty_model
    double E = 1e6;   // Elasticity_modulus
    double nu = 0.495;   // Poisson_ratio
    double pen = 4e9;   // Penalty_parameter
    double C01 = 0.1;   // additional parameter to C10 (optional)
    //weights defined in function - see mat_models_carray.h

    // Step 2: construct test object
    UnitTestIso CANN(matType, E, nu, volType, pen, C01);
    auto &dmn = CANN.com_mod.mockEq.mockDmn;
    auto &stm = dmn.stM;
    auto &w = stm.w;
/*
    //reading from file
    //   FILE* myfile = std::fopen("/Users/divya/svFSIplus/Code/Source/svFSI/weights.txt","r");
	ifstream file("/Users/divya/svFSIplus/Code/Source/svFSI/weights.txt");
    if(!file.is_open()) {
        cerr << "Fail! " << endl;
        // return 1;
    }
    //std::cout << "contents of file:\n" << std::ifstream("/Users/divya/svFSIplus/Code/Source/svFSI/weights.csv").rdbuf() ;

        //Read number using the extraction (>>) operator
        for (int i = 0; i < 2; i++)
        {
          w.push_back({});
          for (int j = 0; j < 16; j++)
          {
            // fscanf(myfile,"%d",&(w[i][j]));
            // std::cout << w[i][j] << std::ends;
            double n=j;
            //file >> n;
            //fscanf(file,"%d",&n);
            std::cout << "n" <<std::endl;
            std::cout << n <<std::endl;
            if (!(file >> n)) {
                cerr << "Fail! 2 " << endl;
                // return 1;
            }
            w[i].push_back(n);
          }
        }
        //Close the file stream
	    //   std::fclose(myfile);
        file.close();*/

    // ifstream file("/Users/divya/svFSIplus/Code/Source/svFSI/ParameterTable.txt");
    // if(!file.is_open()) {
    //     cerr << "Failed to open file!" << endl;
    //     // return 1;
    // }

    // //Read numbers using the extraction (>>) operator
    // for (int i = 0; i < 4; i++) {
    //     w.push_back({});
    //     for (int j = 0; j < 7; j++) {
    //         if (j<=3){
    //             int n;
    //             file >> n;
    //             //cout << n << endl;
    //             w[i].push_back(n);
    //         }
    //         else {
    //             double n;
    //             file >> n;
    //             //cout << n << endl;
    //             w[i].push_back(n);
    //         }
    //         //double n;
    //         //file >> n;
    //         //cout << n << endl;
    //         //w[i].push_back(n);
    //         // cout << j << endl;
    //         // cout << w[i][j] << endl;
    //     }
    // }

    // // Close the file stream
    // file.close();
      


    // Step 3: define the input 
    double F[3][3] = {};
    F[0][0] = 1.0; F[1][1] = 1.0; F[2][2] = 1.0;   // set to Identity

    // Step 4: define the reference output 
    double S_ref[3][3] = {};
    double Dm_ref[6][6] = {};

    // Step 5: run unit test
    CANN.runUnitTest(F, S_ref, Dm_ref);
      
}