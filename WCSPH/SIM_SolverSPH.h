/*
 * Copyright (c) 2024
 *	Side Effects Software Inc.  All rights reserved.
 *
 * Redistribution and use of Houdini Development Kit samples in source and
 * binary forms, with or without modification, are permitted provided that the
 * following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. The name of Side Effects Software may not be used to endorse or
 *    promote products derived from this software without specific prior
 *    written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SIDE EFFECTS SOFTWARE `AS IS' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN
 * NO EVENT SHALL SIDE EFFECTS SOFTWARE BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *----------------------------------------------------------------------------
*/


#ifndef __SIM_SolverSPH_h__
#define __SIM_SolverSPH_h__
#include <UT/UT_Debug.h>

#include <SIM/SIM_SingleSolver.h>
#include <SIM/SIM_OptionsUser.h>
#include <SIM/SIM_Utils.h>


#define SIM_NAME_KERNEL_RADIUS	      "kernel_radius"
#define SIM_NAME_PRESSURE_STRENGTH	  "p_strength"
#define SIM_NAME_GAMMA        	      "gamma"
#define SIM_NAME_REST_DENSITY        	"rest_density"
#define SIM_NAME_DYNAMIC_VISC       	"dynamic_visc"





class SIM_Object;
class SIM_ObjectArray;
class SIM_GeometryCopy;
class SIM_Random;
class GU_NeighbourList;
namespace HDK_Sample {

class SIM_SolverSPH : public SIM_SingleSolver, 
                       public SIM_OptionsUser
{
 public:
   GETSET_DATA_FUNCS_F(SIM_NAME_KERNEL_RADIUS, KernelRadius);
   GETSET_DATA_FUNCS_F(SIM_NAME_PRESSURE_STRENGTH, P_Strength);
   GETSET_DATA_FUNCS_F(SIM_NAME_GAMMA, Gamma);
   GETSET_DATA_FUNCS_F(SIM_NAME_REST_DENSITY, RestDensity);
   GETSET_DATA_FUNCS_F(SIM_NAME_DYNAMIC_VISC, DynamicVisc);


 protected:
   explicit		 SIM_SolverSPH(const SIM_DataFactory *factory); //factory create destory SIM_Data
   ~SIM_SolverSPH() override;
   // This implements your own solver step
   SIM_Result solveSingleObjectSubclass(
               SIM_Engine& engine, SIM_Object& object,
               SIM_ObjectArray& feedback_to_objects,
               const SIM_Time& time_step,
               bool object_is_new) override;

   SIM_Random		*createRandomData(SIM_Object *obj) const; //SIM_Object data for object in sim
   int			 rand_choice(int numchoice, SIM_Random *rand) const;

   float     weight_kernel(const UT_Vector3& a, const UT_Vector3& b);
   UT_Vector3     grad_weight_kernel(const UT_Vector3& a, const UT_Vector3& b);

   void      compute_densities(GU_Detail* gdp, GU_NeighbourList& neighbors);
   void      compute_nonpressure(GU_Detail* gdp, GU_NeighbourList& neighbors);
   void      compute_pressure(GU_Detail* gdp, GU_NeighbourList& neighbors);
   void      clear_accelerations(GU_Detail* gdp);
   void      integrate(GU_Detail* gdp, double dt);
   void      init_attrib(GU_Detail* gdp);

 private:
   static const SIM_DopDescription* getSolverSPHDopDescription();
   float h = 0.1; //default to water
   float rest_density = 1000;
   float p_strength = 50000;
   float gamma = 7;
   float dynamic_viscosity = 0.01;
   DECLARE_STANDARD_GETCASTTOTYPE();
   DECLARE_DATAFACTORY(SIM_SolverSPH,
                      SIM_SingleSolver, 
                       "HDK_SolverSPH",
                       getSolverSPHDopDescription());
};

} // End HDK_Sample namespace

#endif

