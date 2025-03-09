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

 #include "SIM_SolverGravity.h"
 #include <UT/UT_DSOVersion.h>
 #include <UT/UT_Debug.h>
 #include <GU/GU_PrimPoly.h>
 #include <SIM/SIM_Random.h>
 #include <SIM/SIM_RandomTwister.h>
 #include <PRM/PRM_Include.h>
 #include <SIM/SIM_DopDescription.h>
 #include <SIM/SIM_GeometryCopy.h>
 #include <SIM/SIM_DataFilter.h>
 #include <SIM/SIM_Object.h>
 #include <SIM/SIM_ObjectArray.h>
 #include <SIM/SIM_Engine.h>
 #include <SIM/SIM_Force.h>
 #include <cmath>  
 #ifndef M_PI
    #define M_PI 3.14159265358979323846
  #endif
 using namespace HDK_Sample;
 

 
void
initializeSIM(void *)
{
    IMPLEMENT_DATAFACTORY(SIM_SolverGravity);
}
 
 SIM_SolverGravity::SIM_SolverGravity(const SIM_DataFactory *factory)
     : BaseClass(factory),
       SIM_OptionsUser(this)
 {
 }
 
 SIM_SolverGravity::~SIM_SolverGravity()
 {
 }
 
 const SIM_DopDescription *
 SIM_SolverGravity::getSolverGravityDopDescription()
 {
     static PRM_Name	 parm_kc(SIM_NAME_KC, "Kc");
     static PRM_Name	 parm_km(SIM_NAME_KM, "Km");
     static PRM_Name	 parm_kca(SIM_NAME_KCA, "Kca");
     static PRM_Name	 parm_budget(SIM_NAME_BUDGET, "Budget");
     static PRM_Name	 parm_r(SIM_NAME_R, "R");
     static PRM_Name	 parm_ramp(SIM_NAME_RAMP, "R_amp");
     static PRM_Name	 parm_theta(SIM_NAME_THETA, "Theta");
     static PRM_Name	 parm_thetaamp(SIM_NAME_THETAAMP, "Theta_amp");

     static PRM_Template	 theTemplates[] = {
     PRM_Template(PRM_FLT,		1, &parm_kc),
     PRM_Template(PRM_FLT,		1, &parm_km),
     PRM_Template(PRM_FLT,		1, &parm_kca),
     PRM_Template(PRM_FLT,		1, &parm_budget),
     PRM_Template(PRM_FLT,		1, &parm_r),
     PRM_Template(PRM_FLT,		1, &parm_ramp),
     PRM_Template(PRM_ANGLE,		1, &parm_theta),
     PRM_Template(PRM_ANGLE_J,		1, &parm_thetaamp),
     PRM_Template() //<-- what dat guy doing
     };
 
     static SIM_DopDescription	 theDopDescription(true,
                            "hdk_solvergravity",
                            "HDK_SolverGravity",
                            SIM_SOLVER_DATANAME,
                            classname(),
                            theTemplates);
 
     return &theDopDescription;
 }
 
 SIM_Random *
 SIM_SolverGravity::createRandomData(SIM_Object *obj) const
 {
     SIM_Random	*rand = 0;
     
     // Create the random data as subdata attached to the solver. First
     // we look for any existing SIM_Random. If none is found, we create
     // a SIM_RandomTwister.
     rand = SIM_DATA_GET(*obj, "Random", SIM_Random);
     if( !rand )
     rand = SIM_DATA_CREATE(*obj, "Random", SIM_RandomTwister, 0);
 
     return rand;
 }
 int
 SIM_SolverGravity::rand_choice(int numchoice, SIM_Random *rand) const
 {
     int choice = rand->choice(numchoice);
     return choice;
 }
float SIM_SolverGravity::range_limiter(float x_ab_len, float r, float r_amp )
{
    if(x_ab_len <= r) return 1;
    else if (r < x_ab_len && x_ab_len <= r+r_amp) 
    {
        if(r_amp == 0) r_amp = 1e-6;
        return 1 - ((x_ab_len-r)/r_amp);
    }
    else return 0;
}

float SIM_SolverGravity::fov_limiter(float cos_ab, float theta, float theta_amp)
{
    // Convert theta and theta_amp from degrees to radians
    float theta_rad = theta * M_PI / 180.0f;
    float theta_amp_rad = theta_amp * M_PI / 180.0f;

    if(cos_ab >= cos(theta_rad)) return 1;
    else if(cos(theta_rad) > cos_ab && cos_ab > cos(theta_rad + theta_amp_rad)) 
    {
        float denom = (cos(theta_rad) - cos(theta_rad + theta_amp_rad));
        if (denom == 0) denom = 1e-6;
        return 1 - ((cos(theta_rad) - cos_ab) / denom);
    }
    else return 0;
}

SIM_Solver::SIM_Result
SIM_SolverGravity::solveSingleObjectSubclass(SIM_Engine & /*engine*/,
                       SIM_Object &object,
                       SIM_ObjectArray &,
                       const SIM_Time &timestep,
                       bool newobject) //newobject returns true when the node is first connected
{

    // using the solver parameters, e.g., getMyOwnAccuracy()
    float Kc = getKc();
    float Km = getKm();
    float Kca = getKca();
    float r = getR();
    float r_amp = getR_amp();
    float theta = getTheta();
    float theta_amp = getTheta_amp();
    float budget = getBudget();
    const fpreal dt = timestep;
    
    SIM_Result result = SIM_SOLVER_FAIL;

    //Get the object's last state before this time step
    const SIM_Geometry*const geometry = object.getGeometry();
    SIM_Random *rand = createRandomData(&object);

    if( newobject )
    {
        if (!geometry)
        {
            //todo: create geo
            return result = SIM_SOLVER_FAIL;
        }
        // extract simulation state from geometry
        GU_ConstDetailHandle c_gdh = geometry->getGeometry(); 
        const GU_Detail *gdp_in = c_gdh.gdp(); //locking is outdated

        //gu_detail for copy geo
        GU_DetailHandle gdh = c_gdh.getWriteableCopy();// or allocate and set, only for new geo?
        GU_Detail* gdp_out = gdh.gdpNC();

        //get attributes from geo
        //init vels fore new object if doesnt have it
        GA_RWAttributeRef vel_attr = gdp_out->findPointAttribute("v");
        if (!vel_attr.isValid()) {
            //std::cout << "creating vel\n";
            vel_attr = gdp_out->addFloatTuple(GA_ATTRIB_POINT, "v", 3);
            vel_attr.setTypeInfo(GA_TYPE_VECTOR);
            GA_RWHandleV3 vel_h(vel_attr.getAttribute());
            if (!vel_h.isValid()) return SIM_SOLVER_FAIL;
            GA_Offset ptoff;
            GA_Iterator it(gdp_out->getPointRange());
            for (;!it.atEnd(); it.advance())
            {
                ptoff = *it;
                int x = rand_choice(10, rand);
                int y = rand_choice(10, rand);
                int z = rand_choice(10, rand);
                vel_h.set(ptoff, UT_Vector3(x,y,z));
            }
        }
        SIM_GeometryCopy* geometry_copy(
            SIM_DATA_CREATE(
                object, SIM_GEOMETRY_DATANAME, SIM_GeometryCopy,
                SIM_DATA_RETURN_EXISTING | SIM_DATA_ADOPT_EXISTING_ON_DELETE
            )
        );
        
        // store the integrated simulation state in geometry_copy
        geometry_copy->setOwnGeometry(gdh);
        result = SIM_SOLVER_SUCCESS;

    }
    else
    {
        if( geometry )
        {
            // extract simulation state from geometry
            GU_ConstDetailHandle c_gdh = geometry->getGeometry(); 
            const GU_Detail *gdp_in = c_gdh.gdp(); //locking is outdated

            //gu_detail for copy geo
            GU_DetailHandle gdh = c_gdh.getWriteableCopy();// or allocate and set, only for new geo?
            GU_Detail* gdp_out = gdh.gdpNC();

            //get attributes from geo
            GA_RWHandleV3 pos_h(gdp_out->getP());
            GA_RWAttributeRef vel_attr = gdp_out->findPointAttribute("v");
            GA_RWHandleV3 vel_h(vel_attr.getAttribute());
            if (!pos_h.isValid()) return SIM_SOLVER_FAIL;
            if (!vel_h.isValid()) return SIM_SOLVER_FAIL;
            UT_DMatrix4 xform;
            SIMgetGeometryTransform(xform, object);
            // integrate simulation state forward by time step
            GA_Offset ptoff;
            GA_FOR_ALL_PTOFF(gdp_out, ptoff)
            {
                UT_Vector4 pos_local4 = UT_Vector4(pos_h.get(ptoff), 1.0f);
                // Convert to world space
                UT_Vector4 pos_world4 = pos_local4 * xform;
                UT_Vector3 pos_world(pos_world4.x(), pos_world4.y(), pos_world4.z());
                UT_Vector3 vel = vel_h.get(ptoff);

                UT_Vector3 center_acc(0,0,0);
                UT_Vector3 coll_acc(0,0,0);
                UT_Vector3 vel_match_acc(0,0,0);

                GA_Offset ptoffb;
                GA_Iterator it(gdp_out->getPointRange());
                for (;!it.atEnd(); it.advance())
                {
                    ptoffb = *it;
                    if (ptoffb == ptoff) continue; // Skip self
                    UT_Vector4 pos_local4b = UT_Vector4(pos_h.get(ptoffb), 1.0f);
                    // Convert to world space
                    UT_Vector4 pos_world4b = pos_local4b * xform;
                    UT_Vector3 pos_worldb(pos_world4b.x(), pos_world4b.y(), pos_world4b.z());
                    UT_Vector3 velb = vel_h.get(ptoffb);
                    //std::cout << velb.x() << '\n';
                    UT_Vector3 x_ba = (pos_worldb - pos_world);
                    float vel_length = vel.length();
                    if(vel_length <=0 ) vel_length=1e-6;
                    float xba_length = x_ba.length();
                    if(xba_length <=0 ) xba_length=1e-6;
                    fpreal32 cos_ab = x_ba.dot(vel) / (xba_length *vel_length);
                    float xba_dot = x_ba.dot(x_ba);
                    if (xba_dot == 0)  xba_dot = 1e-6;

                    center_acc += x_ba * Kc * range_limiter(x_ba.length(), r, r_amp) * fov_limiter(cos_ab, theta, theta_amp);
                    vel_match_acc += (velb - vel) * Km * range_limiter(x_ba.length(), r, r_amp) * fov_limiter(cos_ab, theta, theta_amp);
                    coll_acc += (x_ba / xba_dot) * -Kca * range_limiter(x_ba.length(), r, r_amp) * fov_limiter(cos_ab, theta, theta_amp);

                }
                UT_Vector3 boid_acc(0,0,0);
                float residual = budget;
                if(residual != 0)
                {
                if( coll_acc.length() < residual ) // Check collision avoidance relative to budget
                {
                    boid_acc += coll_acc;
                    residual -= coll_acc.length();
                    if( vel_match_acc.length() < residual ) // Check velocity matching relative to remaining budget
                    {
                        boid_acc += vel_match_acc;
                        residual -= vel_match_acc.length();
                        if( center_acc.length() < residual ) // Check centering relative to remaining budget
                        {
                            boid_acc += center_acc;
                        }
                        else
                        {
                            boid_acc += center_acc*( residual/center_acc.length() ); // centering exceeds remaining budget, scale it down
                        }
                    }
                    else
                    {
                        boid_acc += vel_match_acc*( residual/vel_match_acc.length() ); // VM exceed remaining budget, scale it down
                    }
                }
                else
                {
                    boid_acc = coll_acc*( residual/coll_acc.length() ); // CA exceeds budget, scale it down
                }
                }

                UT_Vector3 acc(0,0,0);
                acc += boid_acc;
                UT_Vector3 force_world(0, 0, 0);
                acc += force_world;

                // Update the velocity in world space
                vel += acc * dt ;
                
                pos_world += vel * dt;

                UT_Vector3 box_position = {0, 0, 0};  
                UT_Vector3 box_size = {100, 100, 100};  
         // Reflect velocity if it collides with box walls
                if (pos_world.x() > box_position.x() + box_size.x()) { 
                    pos_world.x() = box_position.x() + box_size.x()-1;  
                    vel.x() *= -1;
                }
                else if (pos_world.x() < box_position.x() - box_size.x()) {
                    pos_world.x() = box_position.x() - box_size.x()+1;
                    vel.x() *= -1;
                }

                if (pos_world.y() > box_position.y() + box_size.y()) {
                    pos_world.y() = box_position.y() + box_size.y()-1;
                    vel.y() *= -1;
                }
                else if (pos_world.y() < box_position.y() - box_size.y()) {
                    pos_world.y() = box_position.y() - box_size.y()+1;
                    vel.y() *= -1;
                }

                if (pos_world.z() > box_position.z() + box_size.z()) {
                    pos_world.z() = box_position.z() + box_size.z()-1;
                    vel.z() *= -1;
                }
                else if (pos_world.z() < box_position.z() - box_size.z()) {
                    pos_world.z() = box_position.z() - box_size.z()+1;
                    vel.z() *= -1;
                }

                
                UT_DMatrix4 invXform = xform;
                invXform.invert();
                UT_Vector4 pos_new_local4 = UT_Vector4(pos_world, 1.0f) * invXform;
                UT_Vector3 pos_new_local(pos_new_local4.x(), pos_new_local4.y(), pos_new_local4.z());

    

                pos_h.set(ptoff, pos_new_local);
                vel_h.set(ptoff, vel); //convert vel to local for other nodes??
                // std::cout << pos_world.x() << ' ' << pos_world.y() << ' ' << pos_world.z() << '\n';
                // std::cout << pos_local.x() << ' ' << pos_local.y() << ' ' << pos_local.z() << '\n';
                // std::cout << "Boid Acceleration: " << boid_acc.x() << ", " << boid_acc.y() << ", " << boid_acc.z() << std::endl;
                // std::cout << "Velocity: " << vel.x() << ", " << vel.y() << ", " << vel.z() << std::endl;
                
            }
            SIM_GeometryCopy* geometry_copy(
                SIM_DATA_CREATE(
                    object, SIM_GEOMETRY_DATANAME, SIM_GeometryCopy,
                    SIM_DATA_RETURN_EXISTING | SIM_DATA_ADOPT_EXISTING_ON_DELETE
                )
            );
            
            // store the integrated simulation state in geometry_copy
            geometry_copy->setOwnGeometry(gdh);
            result = SIM_SOLVER_SUCCESS;
        }
    }

    return result;
    
}

 
 