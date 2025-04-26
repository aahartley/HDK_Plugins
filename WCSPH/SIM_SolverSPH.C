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

 #include "SIM_SolverSPH.h"
 #include <UT/UT_DSOVersion.h>
 #include <UT/UT_Debug.h>
 #include <UT/UT_ParallelUtil.h>
 #include <GU/GU_PrimPoly.h>
 #include <GU/GU_NeighbourList.h>
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
 #include <GA/GA_PageIterator.h>
 #include <GA/GA_PageHandle.h>
 using namespace HDK_Sample;
 # define M_PI           3.14159265358979323846  /* pi */


 
void
initializeSIM(void *)
{
    IMPLEMENT_DATAFACTORY(SIM_SolverSPH);
}
 
SIM_SolverSPH::SIM_SolverSPH(const SIM_DataFactory *factory)
     : BaseClass(factory),
       SIM_OptionsUser(this)
 {
 }
 
 SIM_SolverSPH::~SIM_SolverSPH()
 {
 }
 
 const SIM_DopDescription *
 SIM_SolverSPH::getSolverSPHDopDescription()
 {
     static PRM_Name	 parm_KernelRadius(SIM_NAME_KERNEL_RADIUS, "Kernel Radius");
     static PRM_Name	 parm_PressureStrength(SIM_NAME_PRESSURE_STRENGTH, "Pressure Strength");
     static PRM_Name	 parm_Gamma(SIM_NAME_GAMMA, "Gamma");
     static PRM_Name	 parm_RestDensity(SIM_NAME_REST_DENSITY, "Rest Density");
     static PRM_Name	 parm_DynamicVisc(SIM_NAME_DYNAMIC_VISC, "Dynamic Viscosity");
     static PRM_Default defaults[]=
     {
         PRM_Default(0.1f),
         PRM_Default(1000.0f),
         PRM_Default(7.0f),
         PRM_Default(1000.0f),
         PRM_Default(0.01f),
         
     };

     static PRM_Template	 theTemplates[] = {
     PRM_Template(PRM_FLT,		1, &parm_KernelRadius, &defaults[0]),
     PRM_Template(PRM_INT,		1, &parm_PressureStrength, &defaults[0]),
     PRM_Template(PRM_INT,		1, &parm_Gamma, &defaults[0]),
     PRM_Template(PRM_INT,		1, &parm_RestDensity, &defaults[0]),
     PRM_Template(PRM_FLT,		1, &parm_DynamicVisc, &defaults[0]),
     PRM_Template() 
     };
 
     static SIM_DopDescription	 theDopDescription(true,
                            "hdk_solversph",
                            "HDK_SolverSPH",
                            SIM_SOLVER_DATANAME,
                            classname(),
                            theTemplates);
 
     return &theDopDescription;
 }
 
 SIM_Random *
 SIM_SolverSPH::createRandomData(SIM_Object *obj) const
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
 SIM_SolverSPH::rand_choice(int numchoice, SIM_Random *rand) const
 {
     int choice = rand->choice(numchoice);
     return choice;
 }


SIM_Solver::SIM_Result
SIM_SolverSPH::solveSingleObjectSubclass(SIM_Engine & /*engine*/,
                       SIM_Object &object,
                       SIM_ObjectArray &,
                       const SIM_Time &timestep,
                       bool newobject) //newobject returns true when the node is first connected
{

    // using the solver parameters, e.g., getMyOwnAccuracy()
    h = getKernelRadius();
    rest_density = getRestDensity();
    p_strength = getP_Strength();
    gamma = getGamma();
    dynamic_viscosity = getDynamicVisc();
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
        // TODO:attrib map should clean this up?
        //init vels for new object if doesnt have it
        GA_RWHandleV3 vel_h(gdp_out->findPointAttribute("v"));
        if (!vel_h.isValid()) 
        {
            vel_h = GA_RWHandleV3(gdp_out->addFloatTuple(GA_ATTRIB_POINT, "v", 3));
            vel_h->setTypeInfo(GA_TYPE_VECTOR);
            if (!vel_h.isValid()) return SIM_SOLVER_FAIL;
        }
        //acc
        GA_RWHandleV3 acc_h(gdp_out->findPointAttribute("acc"));
        if (!acc_h.isValid()) 
        {
            acc_h = GA_RWHandleV3(gdp_out->addFloatTuple(GA_ATTRIB_POINT, "acc", 3));
            acc_h->setTypeInfo(GA_TYPE_VECTOR);
            if (!acc_h.isValid()) return SIM_SOLVER_FAIL;
        }
        //mass
        GA_RWHandleF mass_h(gdp_out->findPointAttribute("mass"));
        if (!mass_h.isValid()) 
        {
            mass_h = GA_RWHandleF(gdp_out->addFloatTuple(GA_ATTRIB_POINT, "mass", 1));
            if (!mass_h.isValid()) return SIM_SOLVER_FAIL;
        }
        //density
        GA_RWHandleF dens_h(gdp_out->findPointAttribute("dens"));
        if (!dens_h.isValid()) 
        {
            dens_h = GA_RWHandleF(gdp_out->addFloatTuple(GA_ATTRIB_POINT, "dens", 1));
            if (!dens_h.isValid()) return SIM_SOLVER_FAIL;
        }

        init_attrib(gdp_out);
        
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
            GA_RWHandleV3 pos_h(gdp_out->getP()); //assuming the artist didnt change local space
            GA_RWAttributeRef vel_attr = gdp_out->findPointAttribute("v");
            GA_RWHandleV3 vel_h(vel_attr.getAttribute());
            GA_RWAttributeRef acc_attr = gdp_out->findPointAttribute("acc");
            GA_RWHandleV3 acc_h(acc_attr.getAttribute());
            GA_RWAttributeRef mass_attr = gdp_out->findPointAttribute("mass");
            GA_RWHandleF mass_h(mass_attr.getAttribute());
            GA_RWAttributeRef dens_attr = gdp_out->findPointAttribute("dens");
            GA_RWHandleF dens_h(dens_attr.getAttribute());
            if (!pos_h.isValid()) return SIM_SOLVER_FAIL;
            if (!vel_h.isValid()) return SIM_SOLVER_FAIL;
            if (!mass_h.isValid()) return SIM_SOLVER_FAIL;
            if (!dens_h.isValid()) return SIM_SOLVER_FAIL;
            if (!acc_h.isValid()) return SIM_SOLVER_FAIL;
            //init neighborhoods
            GU_NeighbourListParms neigh_parms;
            neigh_parms.setMode(neigh_parms.UNIFORM);
            neigh_parms.setOverrideRadius(true);
            neigh_parms.setRadius(h);
            neigh_parms.setRadiusScale(1);
            GU_NeighbourList neighbors;
            neighbors.build(gdp_out, neigh_parms);
            // SPH solve
            clear_accelerations(gdp_out);
            compute_densities(gdp_out, neighbors);
            compute_nonpressure(gdp_out, neighbors);
            compute_pressure(gdp_out, neighbors);
            integrate(gdp_out, dt);
            neighbors.clearLists(); //needed?
          
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

float SIM_SolverSPH::weight_kernel(const UT_Vector3& a, const UT_Vector3& b)
{
    float weight = 0;
    UT_Vector3 dist = a - b;
    float mag = dist.length();
    float h3 = h * h * h;
    float m_k = static_cast<float>(8.0) / (M_PI*h3);
    const float q = mag/h;
    if (q <= 1.0)
    {
        if (q <= 0.5)
        {
            const float q2 = q*q;
            const float q3 = q2*q;
            weight = m_k * (static_cast<float>(6.0)*q3 - static_cast<float>(6.0)*q2 + static_cast<float>(1.0));
        }
        else
        {
            weight = m_k * (static_cast<float>(2.0)*pow(static_cast<float>(1.0) - q, static_cast<float>(3.0)));
        }
    }
    return weight;
}
UT_Vector3 SIM_SolverSPH::grad_weight_kernel(const UT_Vector3& a, const UT_Vector3& b)
{
    UT_Vector3 weight(0,0,0);
    UT_Vector3 dist = a - b;
    float mag = dist.length();
    float h3 = h * h * h;
    float m_l = static_cast<float>(48.0) / (M_PI*h3);
    const float q = mag/h;
    if ((mag > 1.0e-9) && (q <= 1.0))
    {
        UT_Vector3 gradq = dist / mag;
        gradq /= h;
        if (q <= 0.5)
        {
            weight = m_l*q*((float) 3.0*q - static_cast<float>(2.0))*gradq;
        }
        else
        {
            const float factor = static_cast<float>(1.0) - q;
            weight = m_l*(-factor*factor)*gradq;
        }
    }
    return weight;
}
void SIM_SolverSPH::compute_densities(GU_Detail* gdp, GU_NeighbourList& neighbors)
{
    GA_RWHandleF mass_h(gdp->findPointAttribute("mass"));
    GA_RWHandleF dens_h(gdp->findPointAttribute("dens"));
    GA_RWHandleV3 pos_h(gdp->getP());
    const GA_Range range(gdp->getPointRange());

    GAparallelForEachPage(range, /* shouldthread */ true, [&](GA_PageIterator pit)
    {
        // Create any thread-local data structures etc here.
        GA_RWPageHandleF mass_ph(mass_h.getAttribute());
        GA_RWPageHandleF dens_ph(dens_h.getAttribute());
        GA_ROPageHandleV3 pos_ph(pos_h.getAttribute());
        GAforEachPageBlock(pit, [&](GA_Offset start, GA_Offset end)
        {
            // Perform any per-page setup required.
            mass_ph.setPage(start);
            dens_ph.setPage(start);
            pos_ph.setPage(start);

            for (GA_Offset ptoff = start; ptoff < end; ptoff++)
            {
                GA_IndexArray ptoff_neighbors;
                neighbors.getNeighbours(gdp->pointIndex(ptoff), gdp, ptoff_neighbors);

                float density = 0;
                for (int i = 0; i < ptoff_neighbors.size(); i++)
                {
                    GA_Offset neighbor_ptoff = ptoff_neighbors[i];
                    density += mass_ph.get(neighbor_ptoff) * weight_kernel(pos_ph.get(ptoff),pos_ph.get(neighbor_ptoff));
                }
                dens_ph.set(ptoff, density);
            }
        });
    });
}
// void SIM_SolverSPH::compute_densities(GU_Detail* gdp, GU_NeighbourList& neighbors)
// {
//     GA_Offset ptoff;
//     GA_Iterator it(gdp->getPointRange());
//     GA_RWHandleF mass_h(gdp->findPointAttribute("mass"));
//     GA_RWHandleF dens_h(gdp->findPointAttribute("dens"));
//     GA_RWHandleV3 pos_h(gdp->getP()); //assuming the artist didnt change local space for now
//     int64 N = gdp->getNumPoints();

//     for (;!it.atEnd(); it.advance())
//     {
//         ptoff = *it;
//         GA_IndexArray ptoff_neighbors;
//         neighbors.getNeighbours(gdp->pointIndex(ptoff), gdp, ptoff_neighbors);
  
//         float density = 0;
//         for(int i=0; i <ptoff_neighbors.size(); i++)
//         {
//             GA_Offset ind = (ptoff_neighbors[i]);
//             density += mass_h.get(ptoff) * weight_kernel(pos_h.get(ptoff), pos_h.get(ind));
//         }
//         dens_h.set(ptoff, density);
//     }

// }
void SIM_SolverSPH::compute_nonpressure(GU_Detail* gdp, GU_NeighbourList& neighbors)
{
    GA_RWHandleF mass_h(gdp->findPointAttribute("mass"));
    GA_RWHandleF dens_h(gdp->findPointAttribute("dens"));
    GA_RWHandleV3 pos_h(gdp->getP());
    GA_RWHandleV3 vel_h(gdp->findPointAttribute("v"));
    GA_RWHandleV3 acc_h(gdp->findPointAttribute("acc"));
    const GA_Range range(gdp->getPointRange());
    const UT_Vector3 gravity(0, -9.81, 0);

    GAparallelForEachPage(range, /* shouldthread */ true, [&](GA_PageIterator pit)
    {
        // Create any thread-local data structures etc here.
        GA_ROPageHandleF mass_ph(mass_h.getAttribute());
        GA_ROPageHandleF dens_ph(dens_h.getAttribute());
        GA_ROPageHandleV3 pos_ph(pos_h.getAttribute());
        GA_ROPageHandleV3 vel_ph(vel_h.getAttribute());
        GA_RWPageHandleV3 acc_ph(acc_h.getAttribute());

        GAforEachPageBlock(pit, [&](GA_Offset start, GA_Offset end)
        {
            // Perform any per-page setup required.
            mass_ph.setPage(start);
            dens_ph.setPage(start);
            pos_ph.setPage(start);
            vel_ph.setPage(start);
            acc_ph.setPage(start);

            for (GA_Offset ptoff = start; ptoff < end; ptoff++)
            {
                GA_IndexArray ptoff_neighbors;
                neighbors.getNeighbours(gdp->pointIndex(ptoff), gdp, ptoff_neighbors);

                UT_Vector3 laplacian(0, 0, 0);
                float densitya = dens_ph.get(ptoff);
                float kinematic_viscosity = dynamic_viscosity / densitya;

                // Explicit Viscosity calculation
                for (int i = 0; i < ptoff_neighbors.size(); i++)
                {
                    GA_Offset neighbor_ptoff = ptoff_neighbors[i];
                
                    float densityb = dens_ph.get(neighbor_ptoff);
                    UT_Vector3 v_ij = vel_ph.get(ptoff) - vel_ph.get(neighbor_ptoff);
                    UT_Vector3 x_ij = pos_ph.get(ptoff) - pos_ph.get(neighbor_ptoff);
                    float x_len_sq = x_ij.dot(x_ij);
                    float h_sq = 0.01f * (h * h);  
                    laplacian += ( (mass_ph.get(neighbor_ptoff) / densityb) *
                                (v_ij.dot(x_ij) / (x_len_sq + h_sq)) ) *
                                grad_weight_kernel(pos_ph.get(ptoff), pos_ph.get(neighbor_ptoff));
                }
                laplacian *= 2.0f * (3.0f + 2.0f);  
                UT_Vector3 visc_acc = laplacian * kinematic_viscosity;
                acc_ph.set(ptoff, visc_acc + gravity);
            }
        });
    });
}
// void SIM_SolverSPH::compute_nonpressure(GU_Detail* gdp, GU_NeighbourList& neighbors)
// {
//     GA_Offset ptoff;
//     GA_Iterator it(gdp->getPointRange());
//     GA_RWHandleF mass_h(gdp->findPointAttribute("mass"));
//     GA_RWHandleF dens_h(gdp->findPointAttribute("dens"));
//     GA_RWHandleV3 pos_h(gdp->getP()); //assuming the artist didnt change local space for now (using transform node)
//     GA_RWHandleV3 vel_h(gdp->findPointAttribute("v")); 
//     GA_RWHandleV3 acc_h(gdp->findPointAttribute("acc")); 
//     UT_Vector3 gravity(0,-10,0);
//     float dynamic_viscosity = 0.01;
//     int64 N = gdp->getNumPoints();
//     for (;!it.atEnd(); it.advance())
//     {
//         ptoff = *it;
//         GA_IndexArray ptoff_neighbors;
//         neighbors.getNeighbours(gdp->pointIndex(ptoff), gdp, ptoff_neighbors);
//         UT_Vector3 nonp_acc(0,0,0);
//         UT_Vector3 laplacian(0,0,0);
//         float densitya = dens_h.get(ptoff);
//         float kinematic_viscosity = dynamic_viscosity / densitya;
//         //calc viscosityy
//         for(int i=0; i <ptoff_neighbors.size(); i++)
//         {
//             GA_Offset ind = (ptoff_neighbors[i]);
//             float densityb = dens_h.get(ind);
//             UT_Vector3 v_ij = vel_h.get(ptoff) - vel_h.get(ind);
//             UT_Vector3 x_ij = pos_h.get(ptoff) - pos_h.get(ind);
//             laplacian += ( (mass_h.get(ind) / densityb) * 
//                         ( (v_ij.dot(x_ij)) / ((x_ij.length() * x_ij.length()) + 0.01 * (h*h)) ) ) * grad_weight_kernel(pos_h.get(ptoff), pos_h.get(ind)); 
       
//         }
//         laplacian = (2*(3+2) * laplacian);
//         UT_Vector3 visc_acc = laplacian * kinematic_viscosity;
//         nonp_acc = visc_acc + gravity;
//         acc_h.set(ptoff, nonp_acc);
//     }

// }
void SIM_SolverSPH::compute_pressure(GU_Detail* gdp, GU_NeighbourList& neighbors)
{
    GA_RWHandleF mass_h(gdp->findPointAttribute("mass"));
    GA_RWHandleF dens_h(gdp->findPointAttribute("dens"));
    GA_RWHandleV3 pos_h(gdp->getP());
    GA_RWHandleV3 acc_h(gdp->findPointAttribute("acc"));
    const GA_Range range(gdp->getPointRange());

    GAparallelForEachPage(range, /* shouldthread */ true, [&](GA_PageIterator pit)
    {
        // Create any thread-local data structures etc here.
        GA_ROPageHandleF mass_ph(mass_h.getAttribute());
        GA_ROPageHandleF dens_ph(dens_h.getAttribute());
        GA_ROPageHandleV3 pos_ph(pos_h.getAttribute());
        GA_RWPageHandleV3 acc_ph(acc_h.getAttribute());

        GAforEachPageBlock(pit, [&](GA_Offset start, GA_Offset end)
        {
            // Perform any per-page setup required.
            mass_ph.setPage(start);
            dens_ph.setPage(start);
            pos_ph.setPage(start);
            acc_ph.setPage(start);

            for (GA_Offset ptoff = start; ptoff < end; ptoff++)
            {
                GA_IndexArray ptoff_neighbors;
                neighbors.getNeighbours(gdp->pointIndex(ptoff), gdp, ptoff_neighbors);

                UT_Vector3 p_acc(0, 0, 0);
                const float densitya = std::fmax(dens_ph.get(ptoff), rest_density);

                // Pressure calculation (EOS)
                for (int i = 0; i < ptoff_neighbors.size(); i++)
                {
                    const GA_Offset neighbor_ptoff = ptoff_neighbors[i];
                    
                    const float densityb = std::fmax(dens_ph.get(neighbor_ptoff), rest_density);
                    const float pa = p_strength * (std::pow(densitya/rest_density, gamma) - 1.0f);
                    const float pb = p_strength * (std::pow(densityb/rest_density, gamma) - 1.0f);

                    p_acc += grad_weight_kernel(pos_ph.get(ptoff), pos_ph.get(neighbor_ptoff)) * 
                            mass_ph.get(ptoff) * 
                            ((pa/(densitya*densitya)) + (pb/(densityb*densityb)));
                }
                acc_ph.set(ptoff, acc_ph.get(ptoff) - p_acc);
            }
        });
    });
}
// void SIM_SolverSPH::compute_pressure(GU_Detail* gdp, GU_NeighbourList& neighbors)
// {
//     GA_Offset ptoff;
//     GA_Iterator it(gdp->getPointRange());
//     GA_RWHandleF mass_h(gdp->findPointAttribute("mass"));
//     GA_RWHandleF dens_h(gdp->findPointAttribute("dens"));
//     GA_RWHandleV3 pos_h(gdp->getP()); //assuming the artist didnt change local space for now (using transform node)
//     GA_RWHandleV3 acc_h(gdp->findPointAttribute("acc")); 
//     int64 N = gdp->getNumPoints();

//     for (;!it.atEnd(); it.advance())
//     {
//         ptoff = *it;
//         GA_IndexArray ptoff_neighbors;
//         neighbors.getNeighbours(gdp->pointIndex(ptoff), gdp, ptoff_neighbors);
//         UT_Vector3 p_acc(0,0,0);
//         //calc press EOS
//         for(int i=0; i <ptoff_neighbors.size(); i++)
//         {
//             GA_Offset ind = ptoff_neighbors[i];
//             float densitya = dens_h.get(ptoff);
//             densitya = std::fmax(densitya, rest_density);
//             float densityb = dens_h.get(ind);
//             densityb = std::fmax(densityb, rest_density);
//             float pa = p_strength * (std::pow(densitya/rest_density,gamma) - 1.0);
//             float pb = p_strength * (std::pow(densityb/rest_density,gamma) - 1.0);
//             p_acc += grad_weight_kernel(pos_h.get(ptoff), pos_h.get(ind)) * mass_h.get(ptoff) 
//                      * ((pa/(densitya*densitya)) + (pb/(densityb*densityb))) ;
//         }
//         acc_h.set(ptoff, acc_h.get(ptoff) - p_acc);
//     }

// }
void SIM_SolverSPH::clear_accelerations(GU_Detail* gdp)
{
    GA_RWHandleV3 acc_h(gdp->findPointAttribute("acc"));
    const GA_Range range(gdp->getPointRange());

    GAparallelForEachPage(range, /* shouldthread */ true, [&](GA_PageIterator pit)
    {
        // Create any thread-local data structures etc here.
        GA_RWPageHandleV3 acc_ph(acc_h.getAttribute());

        GAforEachPageBlock(pit, [&](GA_Offset start, GA_Offset end)
        {
            // Perform any per-page setup required.
            acc_ph.setPage(start);

            for (GA_Offset ptoff = start; ptoff < end; ptoff++)
            {
                acc_ph.set(ptoff, UT_Vector3(0,0,0));
            }
        });
    });
}
void SIM_SolverSPH::integrate(GU_Detail* gdp, double dt)
{
    GA_RWHandleV3 pos_h(gdp->getP());
    GA_RWHandleV3 vel_h(gdp->findPointAttribute("v"));
    GA_RWHandleV3 acc_h(gdp->findPointAttribute("acc"));
    const GA_Range range(gdp->getPointRange());

    GAparallelForEachPage(range, /* shouldthread */ true, [&](GA_PageIterator pit)
    {
        // Create any thread-local data structures etc here.
        GA_ROPageHandleV3 pos_ph(pos_h.getAttribute());
        GA_ROPageHandleV3 vel_ph(vel_h.getAttribute());
        GA_RWPageHandleV3 acc_ph(acc_h.getAttribute());

        GAforEachPageBlock(pit, [&](GA_Offset start, GA_Offset end)
        {
            // Perform any per-page setup required.
            acc_ph.setPage(start);
            vel_ph.setPage(start);
            pos_ph.setPage(start);

            for (GA_Offset ptoff = start; ptoff < end; ptoff++)
            {
                UT_Vector3 vel = vel_h.get(ptoff) + (acc_h.get(ptoff)*dt);
                UT_Vector3 pos_world = pos_h.get(ptoff) + (vel_h.get(ptoff)*dt);
                
                vel_ph.set(ptoff, vel);
                pos_ph.set(ptoff, pos_world);            
            }
        });
    });
}
void SIM_SolverSPH::init_attrib(GU_Detail* gdp)
{
    GA_RWHandleV3 pos_h(gdp->getP());
    GA_RWHandleV3 vel_h(gdp->findPointAttribute("v"));
    GA_RWHandleV3 acc_h(gdp->findPointAttribute("acc"));
    GA_RWHandleF dens_h(gdp->findPointAttribute("dens"));
    GA_RWHandleF mass_h(gdp->findPointAttribute("mass"));
    float particle_radius = h/4.0;
    float p_diameter = particle_radius *2;
    const GA_Range range(gdp->getPointRange());
    
    GAparallelForEachPage(range, /* shouldthread */ true, [&](GA_PageIterator pit)
    {
        // Create any thread-local data structures etc here.
        GA_RWPageHandleV3 acc_ph(acc_h.getAttribute());
        GA_RWPageHandleV3 pos_ph(pos_h.getAttribute());
        GA_RWPageHandleV3 vel_ph(vel_h.getAttribute());
        GA_RWPageHandleF dens_ph(dens_h.getAttribute());
        GA_RWPageHandleF mass_ph(mass_h.getAttribute());

        GAforEachPageBlock(pit, [&](GA_Offset start, GA_Offset end)
        {
            // Perform any per-page setup required.
            acc_ph.setPage(start);
            pos_ph.setPage(start);
            vel_ph.setPage(start);
            dens_ph.setPage(start);
            mass_ph.setPage(start);

            for (GA_Offset ptoff = start; ptoff < end; ptoff++)
            {
                // int x = rand_choice(2, rand);
                // int y = rand_choice(5, rand);
                // int z = rand_choice(2, rand);
                acc_ph.set(ptoff, UT_Vector3(0,0,0));
                vel_ph.set(ptoff, UT_Vector3(0,0,0));
                mass_ph.set(ptoff, (p_diameter*p_diameter*p_diameter*0.8f)*rest_density); //volume * rest_dens
                dens_ph.set(ptoff, 0);
            }
        });
    });
}