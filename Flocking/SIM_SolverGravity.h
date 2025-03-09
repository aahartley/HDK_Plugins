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


 #ifndef __SIM_SolverGravity_h__
 #define __SIM_SolverGravity_h__
 
 #include <SIM/SIM_SingleSolver.h>
 #include <SIM/SIM_OptionsUser.h>
 #include <SIM/SIM_Utils.h>



 #define SIM_NAME_KC	      "center_acc_kc"
 #define SIM_NAME_KM	      "vel_match_acc_km"
 #define SIM_NAME_KCA	      "coll_acc_kca"
 #define SIM_NAME_BUDGET	  "budget"
 #define SIM_NAME_R	        "r"
 #define SIM_NAME_RAMP	    "r_amp"
 #define SIM_NAME_THETA	    "theta"
 #define SIM_NAME_THETAAMP	"theta_amp"

 class SIM_Object;
 class SIM_ObjectArray;
 class SIM_GeometryCopy;
 class SIM_Random;

namespace HDK_Sample {
 
class SIM_SolverGravity : public SIM_SingleSolver, 
                        public SIM_OptionsUser
{
  public:
    GETSET_DATA_FUNCS_F(SIM_NAME_KC, Kc);
    GETSET_DATA_FUNCS_F(SIM_NAME_KM, Km);
    GETSET_DATA_FUNCS_F(SIM_NAME_KCA, Kca);
    GETSET_DATA_FUNCS_F(SIM_NAME_BUDGET, Budget);
    GETSET_DATA_FUNCS_F(SIM_NAME_R, R);
    GETSET_DATA_FUNCS_F(SIM_NAME_RAMP, R_amp);
    GETSET_DATA_FUNCS_F(SIM_NAME_THETA, Theta);
    GETSET_DATA_FUNCS_F(SIM_NAME_THETAAMP, Theta_amp);

  protected:
    explicit		 SIM_SolverGravity(const SIM_DataFactory *factory); //factory create destory SIM_Data
    ~SIM_SolverGravity() override;
    // This implements your own solver step
    SIM_Result solveSingleObjectSubclass(
                SIM_Engine& engine, SIM_Object& object,
                SIM_ObjectArray& feedback_to_objects,
                const SIM_Time& time_step,
                bool object_is_new) override;
    float range_limiter(float x_ab_len, float r, float r_amp );
    float fov_limiter(float cos_ab, float theta, float theta_amp);
    SIM_Random		*createRandomData(SIM_Object *obj) const; //SIM_Object data for object in sim
    int			 rand_choice(int numchoice, SIM_Random *rand) const;

  private:
    static const SIM_DopDescription* getSolverGravityDopDescription();
    DECLARE_STANDARD_GETCASTTOTYPE();
    DECLARE_DATAFACTORY(SIM_SolverGravity,
                       SIM_SingleSolver, 
                        "HDK_SolverGravity",
                        getSolverGravityDopDescription());
};
 
 } // End HDK_Sample namespace
 
 #endif
 
 