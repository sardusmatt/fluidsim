/****************************************************************
* FluidSolver - C++ port of Stam's solver						 *
* Author: Matteo Tanca [matteotanca@incomingconnection.net]		  *
*******************************************************************/

#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <iostream>
#include "glew.h"
#include "glut.h" 

#include "MathAux.h"

namespace RTP {

	struct MarkerParticle2D {
		float x;
		float y;
	};

	class FluidSolver {

		protected:
			// Quick matrix implementations
			std::vector<std::vector<float> > _uField;	//u
			std::vector<std::vector<float> > _vField;	//v
			std::vector<std::vector<float> > _prevUField;	//u_prev
			std::vector<std::vector<float> > _prevVField;	//v_prev
			std::vector<std::vector<float> > _densityField;	//dens
			std::vector<std::vector<float> > _prevDensityField;	//dens_prev
			/* Used for vorticity confinement */
			std::vector<std::vector<float> >  _vorticityField;

			bool created;
			bool initialized;

			bool generateRandomVelocities;
			bool confineVorticity;


			/***************************************************************
			****************************************************************
					STABLE FLUIDS VARIABLES + VORTICITY CONFINEMENT
			****************************************************************/
			unsigned int N;
			unsigned int size;
			float dt, diff, visc, vorticityConfinementEps, buoyancyAlphaFactor, buoyancyBetaFactor, ambientTemp;

			float force, source;
			/* Marker particle extension */
			#define NUMOFMARKERPARTICLES 1000
			bool lerpMarkerParticleVelocity;

			MarkerParticle2D particles[NUMOFMARKERPARTICLES];

			// Debug purposes, just to get the center area populated with density
			// when enabled
			bool populateCenter;

		public:

			FluidSolver () {
				created = false;
				initialized = false;
			}

			FluidSolver (int dim, bool activateVorticityConfinement = false, bool populateCenter = false, bool simulateRandomVelocities = false);

			void drawGrid ();

			void drawVelocity ();
			void drawDensity ();

			void draw2DVelocityIn3D ();
			void draw2DDensityIn3D ();

			void clearCurrentMatrices();
			void clearPrevMatrices();

			void setSimulationParameters (float dt = 0.1f, float diff = 0.0f, float visc = 0.0f, float force = 2.0f, float source = 70.0f, float vorticityConfEpsilon = 0.6f, float buoyancyAlphaFactor = 0.0006f, float buoyancyBetaFactor = 0.02f, float ambientTemp = 0.0f);
			
			/* Randomly initializes marker particles, spreading them on the grid */
			void initParticles();

			// Computes the new position of each particle in the marker set
			void updateMarkerParticleSet();

			void drawMarkerParticleSet();

			void switchLERPOnMarker (bool switchOn);

			void clearMatrices();


			void addSimulationInputs (float dx, float dy, float ux, float uy, float vx, float vy, float us, float vs);


			void simulateStep();

			///////////// TEST ///////////////
			void getFromUI (bool mouse_down[3], int win_x, int win_y, int& mx, int& my, int& omx, int& omy);

			void addSource (std::vector<std::vector<float> >& sources, float dt, std::vector<std::vector<float> >& destinationField);


			void setBoundary (int N, int b, std::vector<std::vector<float> >& field);


			void linearSolve (int N, int b, std::vector<std::vector<float> >& field, std::vector<std::vector<float> >& prevField, float a, float c);


			void diffuse (int N, int b, std::vector<std::vector<float> >& field, std::vector<std::vector<float> >& prevField, float diff, float dt);


			void advect (int N, int b, std::vector<std::vector<float> >& dField, std::vector<std::vector<float> >& prevDField,  std::vector<std::vector<float> >& uField, std::vector<std::vector<float> >& vField, float dt);


			void project (int N, std::vector<std::vector<float> >& uField, std::vector<std::vector<float> >& vField, std::vector<std::vector<float> >& pField, std::vector<std::vector<float> >& divField);

			void vorticityConfinement(float dt);

			void addBuoyancy (std::vector<std::vector<float> >& densityField, std::vector<std::vector<float> >& vField, float dt);

			void densityStep ();


			void velocityStep ();


			bool isReady () const {
				return (created && initialized);
			}

			float getDiff () const {
				return diff;
			}

			float getViscosity () const {
				return visc;
			}

			float getTimestep () const {
				return dt;	
			}
			
			float getAmbientTemperature () const {
				return ambientTemp;	
			}
			
			float getBuoyancyAlphaFactor () const {
				return buoyancyAlphaFactor;	
			}
			
			float getBuoyancyBetaFactor () const {
				return buoyancyBetaFactor;
			}
			
			float getVorticityEps () const {
				return vorticityConfinementEps;	
			}

			bool isVorticityConfinementOn () const {
				return confineVorticity;
			}
	};


	////////////////////////////
	/// 3D extension  //////////
	////////////////////////////
	struct MarkerParticle3D {
		float x;
		float y;
		float z;
	};

	/* 3D extension of Stam's solver (FluidSolver) */
	class FluidSolver3D {

		protected:
			// Quick matrix implementations
			std::vector<std::vector<std::vector<float> > >  _uField;	//u
			std::vector<std::vector<std::vector<float> > >  _vField;	//v
			std::vector<std::vector<std::vector<float> > >  _wField;
			std::vector<std::vector<std::vector<float> > >  _prevUField;	//u_prev
			std::vector<std::vector<std::vector<float> > >  _prevVField;	//v_prev
			std::vector<std::vector<std::vector<float> > >  _prevWField;
			std::vector<std::vector<std::vector<float> > >  _densityField;	//dens
			std::vector<std::vector<std::vector<float> > >  _prevDensityField;	//dens_prev
			/* Used for vorticity confinement */
			std::vector<std::vector<std::vector<float> > >  _vorticityField;

			bool created;
			bool initialized;

			bool generateRandomVelocities;
			bool confineVorticity;

			/* In Stam's code is set to 20, as convergence improvements becomes
			 * less and less significant with more iterations, so the computational
			 * effort is not worth the increase in quality of the solution
			 */
			int iterationsInLinearSolver;

			/***************************************************************
			****************************************************************
					STABLE FLUIDS VARIABLES + VORTICITY CONFINEMENT
			****************************************************************/
			unsigned int N;
			unsigned int size;
			float dt, diff, visc, vorticityConfinementEps, buoyancyAlphaFactor, buoyancyBetaFactor, ambientTemp;
			float force, source;

			#define NUMOF3DMARKERPARTICLES 2000
			/* Marker particle extension */
			bool lerpMarkerParticleVelocity;
			MarkerParticle3D particles[NUMOF3DMARKERPARTICLES];

			// Debug purposes, just to get the center area populated with density
			// when enabled
			bool populateCenter;

			RTA::Vector3 origin;

			GLuint textureID;

		public:

			FluidSolver3D () {
				created = false;
				initialized = false;
			}

			FluidSolver3D (int dim, const RTA::Vector3& cubeOrigin = RTA::Vector3(), bool activateVorticityConfinement = false, bool populateCenter = false, bool simulateRandomVelocities = false);
		
			// TODO to be used in conjunction with drawDensity to draw 
			// the smoke with a luminance 3D texture (not implemented yet)
			GLuint get3DTexture ();

			void drawVelocity ();

			void drawContainer ();

			// TODO Assumes drawing is done with a 3dtexture in the pshader */
			void drawDensity ();
			void drawDensityWithBoxes ();

			void clearCurrentMatrices();
			void clearPrevMatrices();

			void setSimulationParameters (float dt = 0.1f, float diff = 0.0f, float visc = 0.0f, float force = 2.0f, float source = 20.0f, float vorticityConfEpsilon = 0.5f, float buoyancyAlphaFactor = 0.0006f, float buoyancyBetaFactor = 0.02f, float ambientTemp = 0.0f, int numOfIterForLinearSolver = 20);

			/* Randomly initializes marker particles, spreading them within the cube */
			void initParticles();

			// Computes the new position of each particle in the marker set
			void updateMarkerParticleSet();

			void drawMarkerParticleSet();

			void switchLERPOnMarker (bool switchOn);

			void clearMatrices();

			void simulateStep();

			///////////// TEST ///////////////
			void getFromUI (bool mouse_down[3], int win_x, int win_y, int& mx, int& my, int& omx, int& omy);

			void addSource (std::vector<std::vector<std::vector<float> > > & sources, float dt, std::vector<std::vector<std::vector<float> > > & destinationField);


			void setBoundary (int N, int b, std::vector<std::vector<std::vector<float> > > & field);


			void linearSolve (int N, int b, std::vector<std::vector<std::vector<float> > > & field, std::vector<std::vector<std::vector<float> > > & prevField, float a, float c);


			void diffuse (int N, int b, std::vector<std::vector<std::vector<float> > > & field, std::vector<std::vector<std::vector<float> > > & prevField, float diff, float dt);


			void advect (int N, int b, std::vector<std::vector<std::vector<float> > > & dField, std::vector<std::vector<std::vector<float> > > & prevDField,  std::vector<std::vector<std::vector<float> > > & uField, std::vector<std::vector<std::vector<float> > > & vField, std::vector<std::vector<std::vector<float> > > & wField, float dt);


			void project (int N, std::vector<std::vector<std::vector<float> > > & uField, std::vector<std::vector<std::vector<float> > > & vField, std::vector<std::vector<std::vector<float> > > & wField, std::vector<std::vector<std::vector<float> > > & pField, std::vector<std::vector<std::vector<float> > > & divField);

			void vorticityConfinement(float dt);

			void addBuoyancy (std::vector<std::vector<std::vector<float> > >& densityField, std::vector<std::vector<std::vector<float> > >& vField, float dt);

			void densityStep ();


			void velocityStep ();


			bool isReady () const {
				return (created && initialized);
			}

			float getDiff () const {
				return diff;
			}

			float getViscosity () const {
				return visc;
			}

			float getTimestep () const {
				return dt;	
			}
			
			float getAmbientTemperature () const {
				return ambientTemp;	
			}
			
			float getBuoyancyAlphaFactor () const {
				return buoyancyAlphaFactor;	
			}
			
			float getBuoyancyBetaFactor () const {
				return buoyancyBetaFactor;
			}
			
			float getVorticityEps () const {
				return vorticityConfinementEps;	
			}

			bool isVorticityConfinementOn () const {
				return confineVorticity;
			}

	};
}

#endif // SOLVER_H

