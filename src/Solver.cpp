#include "Solver.h"

namespace RTP {

	FluidSolver::FluidSolver (int dim, bool activateVorticityConfinement, bool populateCenter, bool simulateRandomVelocities) {
		this->N = dim-2;
		size = (N+2);//*(N+2);
		this->populateCenter = populateCenter;
		clearMatrices();
		
		initParticles();
		created = true;
		confineVorticity = activateVorticityConfinement;
		generateRandomVelocities = simulateRandomVelocities;

	}

	void FluidSolver::clearPrevMatrices () {

		_prevUField.clear();
		_prevUField.resize(size);
		for (unsigned int i = 0; i < size; i++)
			_prevUField[i].resize(size);

		_prevVField.clear();
		_prevVField.resize(size);
		for (unsigned int i = 0; i < size; i++)
			_prevVField[i].resize(size);

		_prevDensityField.clear();
		_prevDensityField.resize(size);
		for (unsigned int i = 0; i < size; i++)
			_prevDensityField[i].resize(size);
	}

	void FluidSolver::clearCurrentMatrices () {
		_uField.clear();
		_uField.resize(size);
		for (unsigned int i = 0; i < size; i++)
			_uField[i].resize(size);

		_vField.clear();
		_vField.resize(size);
		for (unsigned int i = 0; i < size; i++)
			_vField[i].resize(size);

		_densityField.clear();
		_densityField.resize(size);
		for (unsigned int i = 0; i < size; i++)
			_densityField[i].resize(size);

		_vorticityField.clear();
		_vorticityField.resize(size);
		for (unsigned int i = 0; i < size; i++)
			_vorticityField[i].resize(size);

	}

	void FluidSolver::clearMatrices () {
		clearPrevMatrices();
		clearCurrentMatrices();
	}

	void FluidSolver::addSimulationInputs (float dx, float dy, float ux, float uy, float vx, float vy, float us, float vs) {
		clearPrevMatrices();
		int i, j;

		if (dx > 0.0f || dy > 0.0f) {
			i = (int)(dx*N+1);
			j = (int)(dy*N+1);
			if (!(i<1 || i>N || j<1 || j>N)) {
				_prevDensityField[i][j] = source;
			}
		}
		if (vs > 0.0f) {
			i = (int)(vx*N+1);
			j = (int)(vy*N+1);
			if (!(i<1 || i>N || j<1 || j>N)) {
				_prevVField[i][j] = force * vs;
			}
		}
		if (us > 0.0f) {
			i = (int)(ux*N+1);
			j = (int)(uy*N+1);
			if (!(i<1 || i>N || j<1 || j>N)) {
				_prevUField[i][j] = force * us;
			}
		}
	}



	/////////   TEST   //////////////
	void FluidSolver::getFromUI (bool mouse_down[3], int win_x, int win_y, int& mx, int& my, int& omx, int& omy) {
		clearPrevMatrices();

		if (mouse_down[0] || mouse_down[2]) {

			int i = (int)((mx /(float)win_x)*N+1);
			int j = (int)(((win_y-my)/(float)win_y)*N+1);

			if ( i<1 || i>N || j<1 || j>N ) return;

			if ( mouse_down[0] ) {
				_prevUField[i][j] = force * (mx-omx);
				_prevVField[i][j] = force * (omy-my);
			}

			if ( mouse_down[2] ) {
				_prevDensityField[i][j] = source;
			}

			omx = mx;
			omy = my;
		}
	}

	void FluidSolver::initParticles () {
		for (unsigned int i = 0; i < NUMOFMARKERPARTICLES; i++) {
			particles[i].x = RTA::RandomFloat(0.0f, ((float)N));
			particles[i].y = RTA::RandomFloat(0.0f, ((float)N));
		}
	}

	void FluidSolver::switchLERPOnMarker (bool switchOn) {
		lerpMarkerParticleVelocity = switchOn;
	}

	void FluidSolver::setSimulationParameters (float dt, float diff, float visc, float force, float source, float vorticityConfEpsilon, float buoyancyAlphaFactor, float buoyancyBetaFactor, float ambientTemp) {
		this->dt = dt;
		this->diff = diff;
		this->visc = visc;
		this->force = force;
		this->source = source;
		vorticityConfinementEps = vorticityConfEpsilon;
		initialized = true;
		this->buoyancyAlphaFactor = buoyancyAlphaFactor;
		this->buoyancyBetaFactor = buoyancyBetaFactor;
		this->ambientTemp = ambientTemp;
	}

	void FluidSolver::simulateStep () {

		/*if (populateCenter) {
			unsigned int middleX, middleY;
			unsigned int range = N / 10;

			middleX = (int)(_prevDensityField.front().size() / 2.0f);
			middleY = (int)(_prevDensityField.size() / 2.0f);

			for (unsigned int i = middleY - range; i <= middleY + range; i+=2)
				for (unsigned int j = middleX - range; j <= middleX + range; j+=2)
					_prevDensityField[i][j] = source;
		}*/

		if (generateRandomVelocities) {
			static int pass = 0;
			static float multiplier = 0.0f;

			unsigned int middleX, middleY;

			middleX = _prevDensityField.front().size() / 2;
			middleY = _prevDensityField.size() / 2;

			//for (unsigned int k = 0; k < 5; k++) {
				int randomDisplacementX = middleX; //middleX + RTA::RandomInt(-(int)(middleX/2.0), (int)(middleX/2.0-1));
				int randomDisplacementY = middleY; //middleY + RTA::RandomInt(-(int)(middleY/2.0), (int)(middleY/2.0-1));

				float randomScaleFactor1 = 1.0f; //RTA::RandomFloat(-20.0f, 20.0f);
				float randomScaleFactor2 = 1.0f; //RTA::RandomFloat(-20.0f, 20.0f);

				// completely vertical
				//_prevUField[randomDisplacementX][randomDisplacementY] += randomScaleFactor1 * force * multiplier; // (0.5f * pass);
				_prevVField[randomDisplacementX][randomDisplacementY] += randomScaleFactor2 * force * multiplier; // (0.5f * pass);

			//}
			if (pass % 9 == 0)
				_prevDensityField[middleX][middleY] += source * 20.0f; //multiplier;

			pass++;
			if (pass % 2 == 0) {
				if (pass % 10 == 0)
					multiplier = 0.0f;
				multiplier += 0.33f;
			}

		//if (generateRandomVelocities) {
		//	unsigned int middleX, middleY ;
		//	static int pass = 0;

		//	middleX = _prevDensityField.front().size() / 2;
		//	middleY = _prevDensityField.size() / 2;

		//	int rangeX = _prevDensityField.front().size() / 5;
		//	int rangeY = _prevDensityField.size() / 5;

		//	// Add some random forces velocities each step
		//	for (unsigned int k = 0; k < 10; k++) {
		//		int randomDisplacementX = middleX + RTA::RandomInt(-rangeX, rangeX);
		//		int randomDisplacementY = middleY + RTA::RandomInt(-rangeY, rangeY);
		//		float randomScaleFactor1 = RTA::RandomFloat(-7.0f, 7.0f);
		//		float randomScaleFactor2 = RTA::RandomFloat(-7.0f, 7.0f);

		//			
		//		_prevUField[randomDisplacementX][randomDisplacementY] = randomScaleFactor1 * force;
		//		_prevVField[randomDisplacementX][randomDisplacementY] = randomScaleFactor2 * force;
		//	}
		//	pass++;
		//}

		}
		velocityStep();
		densityStep();
		updateMarkerParticleSet();

	}


	// Computes the new position of each particle in the marker set
	void FluidSolver::updateMarkerParticleSet () {
		// Hack to make particle move faster based on the size of the grid
		float velocityScale = logf(N)*logf(N);
		

		for (unsigned int i = 0; i < NUMOFMARKERPARTICLES; i++) {
			
			int cellX = (int)(particles[i].x);
			int cellY = (int)(particles[i].y);

			/* Velocities at the cell */
			float dx = _uField[cellX][cellY];
			float dy = _vField[cellX][cellY];

			if (lerpMarkerParticleVelocity) {
				/* Then determine in which part of the cell the particle is */
				bool firstHalfX = (particles[i].x - cellX) < 0.5f ? true : false;
				bool firstHalfY = (particles[i].y - cellY) < 0.5f ? true : false;

				int v0, h0;
				float vFactor, hFactor;

				if (firstHalfX) {
					v0 = MIN(N, cellX - 1);
					vFactor = -1.0f * (particles[i].x - cellX);
				}
				else {
					v0 = MIN(N, cellX + 1);
					vFactor = 1.0f * (particles[i].x - cellX);
				}

				if (firstHalfY) {
					h0 = MIN(N, cellY - 1);
					hFactor = -1.0f * (particles[i].y - cellY);
				}
				else {
					h0 = MIN(N, cellY + 1);
					hFactor = 1.0f * (particles[i].y - cellY);
				}

				/* Velocities at the neighbourhood */
				float dxv0 = _uField[v0][cellY];
				float dxh0 = _uField[cellX][h0];
				float dxv0h0 = _uField[v0][h0];

				float dyv0 = _vField[v0][cellY];
				float dyh0 = _vField[cellX][h0];
				float dyv0h0 = _vField[v0][h0];

				/* Interpolates velocities */
				dx = LERP(LERP(dx, dxv0, hFactor), LERP(dxh0, dxv0h0, hFactor), vFactor);
				dy = LERP(LERP(dy, dyv0, hFactor), LERP(dyh0, dyv0h0, hFactor), vFactor);

			}
			particles[i].x += dx * velocityScale * dt;
			particles[i].y += dy * velocityScale * dt;

			/* If we are past the boundary reposition randomly */
			if (particles[i].x < 0.0f || particles[i].x > N ||
				particles[i].y < 0.0f || particles[i].y > N) {
				particles[i].x = RTA::RandomFloat(0.0f, ((float)N));
				particles[i].y = RTA::RandomFloat(0.0f, ((float)N));
			}
		}
	}
	

		
	
	
	void FluidSolver::addSource (std::vector<std::vector<float> >& sources, float dt, std::vector<std::vector<float> >& destinationField) {	
		if (sources.size() != destinationField.size() ||
			sources.front().size() != destinationField.front().size()) {
			std::cout << "FluidSolver::addSource - sources.size does not match solver size" << std::endl;
		}
		else {
			for (unsigned int i = 0; i < destinationField.size(); i++)
				for (unsigned j = 0; j < destinationField.front().size(); j++)
					destinationField[i][j] += dt * sources[i][j];
		}
	}

	void FluidSolver::setBoundary (int N, int b, std::vector<std::vector<float> >& field) {
		for (unsigned int i = 1; i <= N; i++) {
			field[0][i] = (b == 1) ? -field[1][i] : field[1][i];
			field[N+1][i] = (b == 1) ? -field[N][i] : field[N][i];
			field[i][0] = (b == 2) ? -field[i][1] : field[i][1];
			field[i][N+1] = (b == 2) ? -field[i][N] : field[i][N];
		}
		field[0][0] = 0.5f * (field[1][0] + field[0][1]);
		field[0][N+1] = 0.5f * (field[1][N+1] + field[0][N]);
		field[N+1][0] = 0.5f * (field[N][0] + field[N+1][1]);
		field[N+1][N+1] = 0.5f * (field[N][N+1] + field[N+1][N]);
	}

	void FluidSolver::linearSolve (int N, int b, std::vector<std::vector<float> >& field, std::vector<std::vector<float> >& prevField, float a, float c) {
		
		for (unsigned int k = 0; k < 20; k++) {
			for (unsigned int i = 1; i <= N; i++) { 
				for (unsigned int j = 1; j <= N; j++) {
					field[i][j] = (prevField[i][j] + a*(field[i-1][j] + field[i+1][j] + field[i][j-1] + field[i][j+1]))/c;

				}
			}
		}
		
		setBoundary (N, b, field);
	}

	void FluidSolver::diffuse (int N, int b, std::vector<std::vector<float> >& field, std::vector<std::vector<float> >& prevField, float diff, float dt) {
		float diffusionrate = dt*diff*N*N;
		linearSolve(N, b, field, prevField, diffusionrate, 1+4*diffusionrate);
	}

	void FluidSolver::advect (int N, int b, std::vector<std::vector<float> >& dField, std::vector<std::vector<float> >& prevDField,  std::vector<std::vector<float> >& uField, std::vector<std::vector<float> >& vField, float dt) {
		int i0, j0, i1, j1;
		float x, y, s0, t0, s1, t1, dt0;

		dt0 = dt*N;

		for (unsigned int i = 1; i <= N; i++) { 
				for (unsigned int j = 1; j <= N; j++) {
					x = i - dt0 * uField[i][j];
					y = j - dt0 * vField[i][j];
					if (x<0.5f) x=0.5f; if (x>N+0.5f) x=N+0.5f; i0=(int)x; i1=i0+1;
					if (y<0.5f) y=0.5f; if (y>N+0.5f) y=N+0.5f; j0=(int)y; j1=j0+1;
					s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;
					dField[i][j] = s0*(t0*prevDField[i0][j0] + t1*prevDField[i0][j1]) + s1*(t0*prevDField[i0][j0] + t1*prevDField[i1][j1]);
				}
		}
		setBoundary(N, b, dField);
	}

	void FluidSolver::project (int N, std::vector<std::vector<float> >& uField, std::vector<std::vector<float> >& vField, std::vector<std::vector<float> >& pField, std::vector<std::vector<float> >& divField) {

		for (unsigned int i = 1; i <= N; i++) { 
				for (unsigned int j = 1; j <= N; j++) {
					divField[i][j] = -0.5f * (uField[i+1][j] - uField[i-1][j] + vField[i][j+1] - vField[i][j-1]) / N;
					pField[i][j] = 0.0f;
				}
		}
			
		setBoundary(N, 0, divField);
		setBoundary(N, 0, pField);

		linearSolve(N, 0, pField, divField, 1, 4);

		for (unsigned int i = 1; i <= N; i++) { 
				for (unsigned int j = 1; j <= N; j++) {
					uField[i][j] -= 0.5f * N * (pField[i+1][j] - pField[i-1][j]);
					vField[i][j] -= 0.5f * N * (pField[i][j+1] - pField[i][j-1]);
				}
		}

		setBoundary(N, 1, uField);
		setBoundary(N, 2, vField);
	}

	/* The vorticity confinement for each cell is computed to preserve turbulences that
	 * are suppressed due to artificial numerical dissipation.
	 * The forces at points (i, j) are given by N x omega, where N are normalized
	 * vorticity location vectors pointing from lower vorticity concentrations to
	 * higher vorticity concentrations and omega is the vorticity.
	 * The forces are the scaled by multiplying by epsilon, which gives control
	 * on amount of small scale detail to be added back in the flow
	 */
	void FluidSolver::vorticityConfinement (float dt) {

		// Uses the prev velocity fields as temporary buffers, since
		// the computation happens before the swap
		
		float dt0 = dt * vorticityConfinementEps;

		float x, y;
		
		for (unsigned int i = 1; i < N; i++) { 
				for (unsigned int j = 1; j < N; j++) {

				x = _prevUField[i][j] = (_vField[i+1][j] - _vField[i-1][j]) * 0.5f;
				y = _prevVField[i][j] = (_uField[i][j+1] - _uField[i][j-1]) * 0.5f;

				_vorticityField[i][j] = sqrtf(x*x+y*y);
			}
		}

		for (unsigned int i = 1; i < N; i++) { 
			for (unsigned int j = 1; j < N; j++) {
				float Nx = (_vorticityField[i+1][j] - _vorticityField[i-1][j]) * 0.5f; 
				float Ny = (_vorticityField[i][j+1] - _vorticityField[i][j-1]) * 0.5f;

				// add small epsilon to avoid / 0
				float lenght = 1.0f/(sqrtf(Nx*Nx+Ny*Ny)+0.0000001f);
				Nx *= lenght;
				Ny *= lenght;

				/* Uses values from the temporary fields */
				_uField[i][j] += Ny * _prevVField[i][j] * dt0;
				_vField[i][j] += Nx * _prevUField[i][j] * dt0;
			}
		}

	}

	/* Implements the buoyancy force computation, as described in "Visual simulation
	 * of smoke" by Stam et al., with two additional assumptions.
	 * 1) Temperature is assumed being a fraction of density at each point.
	 * (This is in general not true. In fact, assuming an ideal gas at constant pressure
	 * the temperature would be inversely proportional to density). Hence the
	 * implementation is simplified a bit, with no need to keep and
	 * update an additional field to track temperature at each voxel.
	 * 2) Ambient temperature is assumed constant over all grid.
	 * The formula to compute the force is:
	 * -alpha rho[i,j] + beta (temp[i,j] - ambientTemperature)
	 * The result is then applied to the v velocity field simply multiplying the
	 * force by timestep (that is, dt parameter), as described in the paper.
	 */
	void FluidSolver::addBuoyancy (std::vector<std::vector<float> >& densityField, std::vector<std::vector<float> >& vField, float dt) {

		for (unsigned int i = 0; i < densityField.size(); i++)
			for (unsigned int j = 0; j < densityField.front().size(); j++) {
				float temp = 0.6 * densityField[i][j];
				vField[i][j] += ((-buoyancyAlphaFactor * densityField[i][j]) + buoyancyBetaFactor * (temp - ambientTemp)) * dt;
			}

	}

	void FluidSolver::densityStep () {
		addSource(_prevDensityField, dt, _densityField);
		_prevDensityField.swap(_densityField);
		diffuse(N, 0, _densityField, _prevDensityField, diff, dt);
		_prevDensityField.swap(_densityField);
		advect(N, 0, _densityField, _prevDensityField, _uField, _vField, dt);
	}

	void FluidSolver::velocityStep () {

		addSource(_prevUField, dt, _uField); 
		addSource(_prevVField, dt, _vField);
	
		/**********************************
		* VORTICITY CONFINEMENT EXTENSION *
		**********************************/
		if (confineVorticity) {
			addBuoyancy(_densityField, _vField, dt);
			vorticityConfinement(dt);
		}

		_prevUField.swap(_uField);
		diffuse(N, 1, _uField, _prevUField, visc, dt);

		_prevVField.swap(_vField);
		diffuse(N, 2, _vField, _prevVField, visc, dt);
		project(N, _uField, _vField, _prevUField, _prevVField);

		_prevUField.swap(_uField);
		_prevVField.swap(_vField);
		advect(N, 1, _uField, _prevUField, _prevUField, _prevVField, dt);
		advect(N, 2, _vField, _prevVField, _prevUField, _prevVField, dt);
		project(N, _uField, _vField, _prevUField, _prevVField);
	}


	void FluidSolver::drawGrid () {

		float h = 1.0f / N;

		glColor3f(0.0f, 1.0f, 1.0f);

		glLineWidth(1.0f);
		glBegin(GL_LINES);
			for (unsigned int i = 0; i <= N; i++) {
				glVertex2f(0.0f, i * h);
				glVertex2f(1.0f, i * h);
				glVertex2f(i * h, 0.0f);
				glVertex2f(i * h, 1.0f);
			}

		glEnd ();

		glColor3f(1.0f, 1.0f, 1.0f);
	}

	void FluidSolver::draw2DDensityIn3D () {
		float x, y, d00, d01, d10, d11;

		float z = 5.0f;

		float h = 2.0f/N;

		glBegin(GL_QUADS);
			for (unsigned int i = 0; i <= N; i++) {
				x = (i - 0.5f) * h;
				for (unsigned int j = 0; j <= N; j++) {
					y = (j - 0.5f) * h;

					d00 = _densityField[i][j];
					d01 = _densityField[i][j+1];
					d10 = _densityField[i+1][j];
					d11 = _densityField[i+1][j+1];

					glColor3f(d00, d00, d00); 
					glVertex3f(x, y, z);

					glColor3f(d10, d10, d10);
					glVertex3f(x+h, y, z);
					glColor3f(d11, d11, d11);
					glVertex3f(x+h, y+h, z);
					glColor3f(d01, d01, d01);
					glVertex3f(x, y+h, z);
				}
			}

		glEnd ();

		glColor3f(1.0f, 1.0f, 1.0f);
	}

	void FluidSolver::drawDensity () {
		float x, y, d00, d01, d10, d11;

		float h = 1.0f / N;

		glBegin(GL_QUADS);
			for (unsigned int i = 0; i <= N; i++) {
				x = (i - 0.5f) * h;
				for (unsigned int j = 0; j <= N; j++) {
					y = (j - 0.5f) * h;

					d00 = _densityField[i][j];
					d01 = _densityField[i][j+1];
					d10 = _densityField[i+1][j];
					d11 = _densityField[i+1][j+1];

					glColor3f(d00, d00, d00); 
					glVertex2f(x, y);

					glColor3f(d10, d10, d10);
					glVertex2f(x+h, y);
					glColor3f(d11, d11, d11);
					glVertex2f(x+h, y+h);
					glColor3f(d01, d01, d01);
					glVertex2f(x, y+h);
				}
			}

		glEnd ();

		glColor3f(1.0f, 1.0f, 1.0f);
	}

	// TODO implement
	void FluidSolver::draw2DVelocityIn3D () {
		float x, y;
		float z = 5.0f;
		float h = 2.0f/N;

		glColor3f(1.0f, 1.0f, 0.0f);
		glLineWidth(1.0f);

		glBegin (GL_LINES);

			for (unsigned int i = 0; i <= N; i++) {
				x = (i - 0.5f) * h;
				for (unsigned int j = 0; j <= N; j++) {
					y = (j - 0.5f) * h;
					glVertex3f(x, y, z);
					glVertex3f(x+_uField[i][j], y+_vField[i][j], z);
				}
			}

		glEnd ();

		glColor3f(1.0f, 1.0f, 1.0f);
	}

	void FluidSolver::drawMarkerParticleSet () {
		
		float h = 1.0f/N;

		glColor3f(1.0f, 0.0f, 0.0f);
		glPointSize(1.0f);

		glBegin (GL_POINTS);

			for (unsigned int i = 0; i < NUMOFMARKERPARTICLES; i++) {

				glVertex2f(particles[i].x * h, particles[i].y * h);
			}

		glEnd ();

		glColor3f(1.0f, 1.0f, 1.0f);
		glPointSize(1.0f);
	}

	void FluidSolver::drawVelocity () {
		float x, y;
		float h = 1.0f/N;

		glColor3f(1.0f, 1.0f, 0.0f);
		glLineWidth(1.0f);

		glBegin (GL_LINES);

			for (unsigned int i = 0; i <= N; i++) {
				x = (i - 0.5f) * h;
				for (unsigned int j = 0; j <= N; j++) {
					y = (j - 0.5f) * h;
					glVertex2f(x, y);
					glVertex2f(x+_uField[i][j], y+_vField[i][j]);
				}
			}

		glEnd ();

		glColor3f(1.0f, 1.0f, 1.0f);


	}

	
	////////////////////////////////////////////
	///////////////		3D extension	////////
	////////////////////////////////////////////
	FluidSolver3D::FluidSolver3D (int dim, const RTA::Vector3& cubeOrigin, bool activateVorticityConfinement, bool populateCenter, bool simulateRandomVelocities) {
		N = dim-2;
		size = (N+2);//*(N+2);
		this->populateCenter = populateCenter;
		clearMatrices();
		created = true;
		origin = RTA::Vector3(cubeOrigin.x, cubeOrigin.y, cubeOrigin.z);
		confineVorticity = activateVorticityConfinement;
		generateRandomVelocities = simulateRandomVelocities;
		initParticles();	/* Initialize marker particle set */
	}

	void FluidSolver3D::clearPrevMatrices () {

		_prevUField.clear();
		_prevUField.resize(size);
		for (unsigned int i = 0; i < size; i++) {
			_prevUField[i].resize(size);
			for (unsigned int j = 0; j < size; j++)
				_prevUField[i][j].resize(size);
		}
				

		_prevVField.clear();
		_prevVField.resize(size);
		for (unsigned int i = 0; i < size; i++) {
			_prevVField[i].resize(size);
			for (unsigned int j = 0; j < size; j++)
				_prevVField[i][j].resize(size);
		}

		_prevWField.clear();
		_prevWField.resize(size);
		for (unsigned int i = 0; i < size; i++) {
			_prevWField[i].resize(size);
			for (unsigned int j = 0; j < size; j++)
				_prevWField[i][j].resize(size);
		}

		_prevDensityField.clear();
		_prevDensityField.resize(size);
		for (unsigned int i = 0; i < size; i++) {
			_prevDensityField[i].resize(size);
			for (unsigned int j = 0; j < size; j++)
				_prevDensityField[i][j].resize(size);
		}
	}

	void FluidSolver3D::clearCurrentMatrices () {
		_uField.clear();
		_uField.resize(size);
		for (unsigned int i = 0; i < size; i++) {
			_uField[i].resize(size);
			for (unsigned int j = 0; j < size; j++)
				_uField[i][j].resize(size);
		}

		_vField.clear();
		_vField.resize(size);
		for (unsigned int i = 0; i < size; i++) {
			_vField[i].resize(size);
			for (unsigned int j = 0; j < size; j++)
				_vField[i][j].resize(size);
		}

		_wField.clear();
		_wField.resize(size);
		for (unsigned int i = 0; i < size; i++) {
			_wField[i].resize(size);
			for (unsigned int j = 0; j < size; j++)
				_wField[i][j].resize(size);
		}

		_densityField.clear();
		_densityField.resize(size);
		for (unsigned int i = 0; i < size; i++) {
			_densityField[i].resize(size);
			for (unsigned int j = 0; j < size; j++)
				_densityField[i][j].resize(size);
		}

		_vorticityField.clear();
		_vorticityField.resize(size);
		for (unsigned int i = 0; i < size; i++) {
			_vorticityField[i].resize(size);
			for (unsigned int j = 0; j < size; j++)
				_vorticityField[i][j].resize(size);
		}
	}

	void FluidSolver3D::clearMatrices () {
		clearPrevMatrices();
		clearCurrentMatrices();
	}

	/////////   TEST   //////////////
	void FluidSolver3D::getFromUI (bool mouse_down[3], int win_x, int win_y, int& mx, int& my, int& omx, int& omy) {
		clearPrevMatrices();

		if (mouse_down[0] || mouse_down[2]) {

			int i = (int)((mx /(float)win_x)*N+1);
			int j = (int)(((win_y-my)/(float)win_y)*N+1);
			int l = (int)(((win_x-mx)/(float)win_x)*N+1);

			if ( i<1 || i>N || j<1 || j>N || l<1 || l>N) return;

			if ( mouse_down[0] ) {
				_prevUField[i][j][l] = force * (mx-omx);
				_prevVField[i][j][l] = force * (omy-my);
				_prevWField[i][j][l] = force * (omy-my + mx-omx)/2.0f;	// TEST
			}

			if ( mouse_down[2] ) {
				_prevDensityField[i][j][l] = source;
			}

			omx = mx;
			omy = my;
		}
	}


	void FluidSolver3D::setSimulationParameters (float dt, float diff, float visc, float force, float source, float vorticityConfEpsilon, float buoyancyAlphaFactor, float buoyancyBetaFactor, float ambientTemp, int numOfIterForLinearSolver) {
		this->dt = dt;
		this->diff = diff;
		this->visc = visc;
		this->force = force;
		this->source = source;
		
		iterationsInLinearSolver = numOfIterForLinearSolver;
		initialized = true;
		vorticityConfinementEps = vorticityConfEpsilon;
		this->buoyancyAlphaFactor = buoyancyAlphaFactor;
		this->buoyancyBetaFactor = buoyancyBetaFactor;
		this->ambientTemp = ambientTemp;
	}

	void FluidSolver3D::switchLERPOnMarker (bool switchOn) {
		lerpMarkerParticleVelocity = switchOn;
	}

	void FluidSolver3D::initParticles () {
		for (unsigned int i = 0; i < NUMOF3DMARKERPARTICLES; i++) {
			particles[i].x = RTA::RandomFloat(0.0f, ((float)N));
			particles[i].y = RTA::RandomFloat(0.0f, ((float)N));
			particles[i].z = RTA::RandomFloat(0.0f, ((float)N));
		}
	}

	void FluidSolver3D::updateMarkerParticleSet () {
		// Hack to make particle move faster
		float velocityScale = pow(logf(N), 3);

		for (unsigned int i = 0; i < NUMOF3DMARKERPARTICLES; i++) {
			
			int cellX = (int)(particles[i].x);
			int cellY = (int)(particles[i].y);
			int cellZ = (int)(particles[i].z);

			/* Velocities at the cell */
			float dx = _uField[cellX][cellY][cellZ];
			float dy = _vField[cellX][cellY][cellZ];
			float dz = _wField[cellX][cellY][cellZ];

			if (lerpMarkerParticleVelocity) {
				/* Then determine in which part of the cell the particle is */
				bool firstHalfX = (particles[i].x - cellX) < 0.5f ? true : false;
				bool firstHalfY = (particles[i].y - cellY) < 0.5f ? true : false;
				bool firstHalfZ = (particles[i].z - cellZ) < 0.5f ? true : false;

				int v0, h0, w0;
				float vFactor, hFactor, wFactor;

				if (firstHalfX) {
					v0 = MIN(N, cellX - 1);
					vFactor = -1.0f * (particles[i].x - cellX);
				}
				else {
					v0 = MIN(N, cellX + 1);
					vFactor = 1.0f * (particles[i].x - cellX);
				}

				if (firstHalfY) {
					h0 = MIN(N, cellY - 1);
					hFactor = -1.0f * (particles[i].y - cellY);
				}
				else {
					h0 = MIN(N, cellY + 1);
					hFactor = 1.0f * (particles[i].y - cellY);
				}

				if (firstHalfZ) {
					w0 = MIN(N, cellZ - 1);
					wFactor = -1.0f * (particles[i].z - cellZ);
				}
				else {
					w0 = MIN(N, cellZ + 1);
					wFactor = 1.0f * (particles[i].z - cellZ);
				}

				/* Velocities at the neighbourhood */
				float dxv0 = _uField[v0][cellY][cellZ];
				float dxh0 = _uField[cellX][h0][cellZ];
				float dxw0 = _uField[cellX][cellY][w0];
				float dxv0h0 = _uField[v0][h0][cellZ];
				float dxh0w0 = _uField[cellX][h0][w0];
				float dxv0w0 = _uField[v0][cellY][w0];
				float dxv0h0w0 = _uField[v0][h0][w0];

				float dyv0 = _vField[v0][cellY][cellZ];
				float dyh0 = _vField[cellX][h0][cellZ];
				float dyw0 = _vField[cellX][cellY][w0];
				float dyv0h0 = _vField[v0][h0][cellZ];
				float dyh0w0 = _vField[cellX][h0][w0];
				float dyv0w0 = _vField[v0][cellY][w0];
				float dyv0h0w0 = _vField[v0][h0][w0];

				float dzv0 = _wField[v0][cellY][cellZ];
				float dzh0 = _wField[cellX][h0][cellZ];
				float dzw0 = _wField[cellX][cellY][w0];
				float dzv0h0 = _wField[v0][h0][cellZ];
				float dzh0w0 = _wField[cellX][h0][w0];
				float dzv0w0 = _wField[v0][cellY][w0];
				float dzv0h0w0 = _wField[v0][h0][w0];

				/* Interpolates velocities */
				dx = LERP(
					(LERP(LERP(dx, dxv0, wFactor), LERP(dxh0, dxw0, vFactor), hFactor)),
					(LERP(LERP(dxv0h0, dxh0w0, wFactor), LERP(dxv0w0, dxv0h0w0, hFactor), vFactor)),
					(LERP(hFactor, wFactor, vFactor)));

				dy = LERP(
					(LERP(LERP(dy, dyv0, wFactor), LERP(dyh0, dyw0, vFactor), hFactor)),
					(LERP(LERP(dyv0h0, dyh0w0, wFactor), LERP(dyv0w0, dyv0h0w0, hFactor), vFactor)),
					(LERP(wFactor, vFactor, hFactor)));

				dz = LERP(
					(LERP(LERP(dz, dzv0, wFactor), LERP(dzh0, dzw0, vFactor), hFactor)),
					(LERP(LERP(dzv0h0, dzh0w0, wFactor), LERP(dzv0w0, dzv0h0w0, hFactor), vFactor)),
					(LERP(vFactor, hFactor, wFactor)));
			}

			particles[i].x += dx * velocityScale * dt;
			particles[i].y += dy * velocityScale * dt;
			particles[i].z += dz * velocityScale * dt;

			/* If we are past the boundary reposition randomly */
			if (particles[i].x < 0.0f || particles[i].x > N ||
				particles[i].y < 0.0f || particles[i].y > N ||
				particles[i].z < 0.0f || particles[i].z > N) {
				particles[i].x = RTA::RandomFloat(0.0f, ((float)N));
				particles[i].y = RTA::RandomFloat(0.0f, ((float)N));
				particles[i].z = RTA::RandomFloat(0.0f, ((float)N));
			}
			
		}
	}

	void FluidSolver3D::simulateStep () {

		/*if (populateCenter) {
			unsigned int middleX, middleY, middleZ;
			unsigned int range = N / 30.0f;

			middleX = (int)(_prevDensityField.front().size() / 2.0f);
			middleY = (int)(_prevDensityField.size() / 2.0f);
			middleZ = (int)(_prevDensityField.front().front().size() / 2.0f);

			for (unsigned int i = middleY - range; i <= middleY + range; i++)
				for (unsigned int j = middleX - range; j <= middleX + range; j++)
					for (unsigned int l = middleZ - range; l <= middleZ + range; l++)
						_prevDensityField[i][j][l] = source;
			
			populateCenter = false;
		}*/

		if (generateRandomVelocities) {
			static int pass = 0;
			static float multiplier = 0.0f;

			unsigned int middleX, middleY, middleZ;

			middleX = _prevDensityField.front().size() / 2;
			middleY = _prevDensityField.size() / 2;
			middleZ = _prevDensityField.front().front().size() / 2;

			//for (unsigned int k = 0; k < 5; k++) {
				int randomDisplacementX = middleX; //middleX + RTA::RandomInt(-(int)(middleX/2.0), (int)(middleX/2.0-1));
				int randomDisplacementY = middleY; //middleY + RTA::RandomInt(-(int)(middleY/2.0), (int)(middleY/2.0-1));
				int randomDisplacementZ = middleZ; //middleZ + RTA::RandomInt(-(int)(middleZ/2.0), (int)(middleZ/2.0-1));
				float randomScaleFactor1 = 1.0f; //RTA::RandomFloat(-20.0f, 20.0f);
				float randomScaleFactor2 = 1.0f; //RTA::RandomFloat(-20.0f, 20.0f);
				float randomScaleFactor3 = 1.0f; //RTA::RandomFloat(-20.0f, 20.0f);

				// only vertical component
				//_prevUField[randomDisplacementX][randomDisplacementY][randomDisplacementZ] += randomScaleFactor1 * force * multiplier; // (0.5f * pass);
				_prevVField[randomDisplacementX][randomDisplacementY][randomDisplacementZ] += randomScaleFactor2 * force * multiplier; // (0.5f * pass);
				//_prevWField[randomDisplacementX][randomDisplacementY][randomDisplacementZ] += randomScaleFactor3 * force * multiplier; // (0.5f * pass);
			//}

			if (pass % 10 == 0)
				//Bottom
				_prevDensityField[middleX][middleY][middleZ] += source * 10.0f; //multiplier;

			pass++;
			if (pass % 2 == 0) {
				multiplier += 0.33f;
				//TEST
			}


		}

		velocityStep();
		densityStep();
		updateMarkerParticleSet();
	}
	
	void FluidSolver3D::addSource (std::vector<std::vector<std::vector<float> > > & sources, float dt, std::vector<std::vector<std::vector<float> > > & destinationField) {	
		if (sources.size() != destinationField.size() ||
			sources.front().size() != destinationField.front().size() ||
			sources.front().front().size() != destinationField.front().front().size()) {
			std::cout << "FluidSolver3D::addSource - sources.size does not match solver size" << std::endl;
		}
		else {
			for (unsigned int i = 0; i < destinationField.size(); i++)
				for (unsigned j = 0; j < destinationField.front().size(); j++)
					for (unsigned l = 0; l < destinationField.front().front().size(); l++)
						destinationField[i][j][l] += dt * sources[i][j][l];
		}
	}

	void FluidSolver3D::setBoundary (int N, int b, std::vector<std::vector<std::vector<float> > > & field) {
		for (unsigned int i = 1; i <= N; i++) {
			for (unsigned int j = 1; j <= N; j++) {
				field[0][i][j] = (b == 1) ? -field[1][i][j] : field[1][i][j];
				field[N+1][i][j] = (b == 1) ? -field[N][i][j] : field[N][i][j];
				field[i][0][j] = (b == 2) ? -field[i][1][j] : field[i][1][j];
				field[i][N+1][j] = (b == 2) ? -field[i][N][j] : field[i][N][j];
				field[i][j][0] = (b == 3) ? -field[i][j][1] : field[i][j][1];
				field[i][j][N+1] = (b == 3) ? -field[i][j][N] : field[i][j][N];
			}
		}
		field[0][0][0] = (field[1][0][0] + field[0][1][0] + field[0][0][1]) / 3.0f;
		field[0][N+1][0] = (field[1][N+1][0] + field[0][N][0] + field[0][N+1][1]) / 3.0f;
		field[0][0][N+1] = (field[1][0][N+1] + field[0][1][N+1] + field[0][0][N]) / 3.0f;
		field[N+1][0][0] = (field[N][0][0] + field[N+1][1][0] + field[N+1][0][1]) / 3.0f;
		field[N+1][N+1][0] = (field[N][N+1][0] + field[N+1][N][0] + field[N+1][N+1][1]) / 3.0f;
		field[0][N+1][N+1] = (field[1][N+1][N+1] + field[0][N][N+1] + field[0][N+1][N]) / 3.0f;
		field[N+1][0][N+1] = (field[N][0][N+1] + field[N+1][1][N+1] + field[N+1][0][N]) / 3.0f;
		field[N+1][N+1][N+1] = (field[N][N+1][N+1] + field[N+1][N][N+1] + field[N+1][N+1][N]) / 3.0f;
	}

	void FluidSolver3D::linearSolve (int N, int b, std::vector<std::vector<std::vector<float> > > & field, std::vector<std::vector<std::vector<float> > > & prevField, float a, float c) {
		
		for (unsigned int k = 0; k < iterationsInLinearSolver; k++) {
			for (unsigned int i = 1; i <= N; i++) { 
				for (unsigned int j = 1; j <= N; j++) {
					for (unsigned int l = 1; l <= N; l++) {
						field[i][j][l] = (prevField[i][j][l] + a*(field[i-1][j][l] + field[i+1][j][l] + field[i][j-1][l] + field[i][j+1][l] + field[i][j][l-1] + field[i][j][l+1]))/c;
					}
					
				}
			}
		}
		
		setBoundary (N, b, field);
	}

	void FluidSolver3D::diffuse (int N, int b, std::vector<std::vector<std::vector<float> > > & field, std::vector<std::vector<std::vector<float> > > & prevField, float diff, float dt) {
		float diffusionrate = dt*diff*N*N*N;
		linearSolve(N, b, field, prevField, diffusionrate, 1+6*diffusionrate);
	}

	void FluidSolver3D::advect (int N, int b, std::vector<std::vector<std::vector<float> > > & dField, std::vector<std::vector<std::vector<float> > > & prevDField,  std::vector<std::vector<std::vector<float> > > & uField, std::vector<std::vector<std::vector<float> > > & vField, std::vector<std::vector<std::vector<float> > > & wField, float dt) {
		int i0, j0, l0, i1, j1, l1;
		float x, y, z, s0, t0, u0, s1, t1, u1, dt0;

		dt0 = dt*N;

		for (unsigned int i = 1; i <= N; i++) { 
			for (unsigned int j = 1; j <= N; j++) {
				for (unsigned int l = 1; l <= N; l++) {
					x = ((float)i) - dt0 * uField[i][j][l];
					y = ((float)j) - dt0 * vField[i][j][l];
					z = ((float)l) - dt0 * wField[i][j][l];

					if (x<0.5f)
						x=0.5f;
					if (x>((float)N)+0.5f) 
						x=((float)N)+0.5f;
					i0=(int)x;
					i1=i0+1;

					if (y<0.5f)
						y=0.5f;
					if (y>((float)N)+0.5f)
						y=((float)N)+0.5f;
					j0=(int)y;
					j1=j0+1;

					if (z<0.5f)
						z=0.5f;
					if (z>((float)N)+0.5f)
						z=((float)N)+0.5f;
					l0=(int)z;
					l1=l0+1;

					s1 = x-((float)i0);
					s0 = 1.0f-s1;

					t1 = y-((float)j0);
					t0 = 1.0f-t1;
					u1 = z-((float)l0);
					u0 = 1.0f-u1;


					// NOT SO SURE THIS IS CORRECT
					/*float temp0 = s0*(t0*prevDField[i0][j0][l0] + t1*prevDField[i0][j1][l0]) + s1*(t0*prevDField[i1][j0][l0] + t1*prevDField[i1][j1][l0]);

					float temp1 = s0*(t0*prevDField[i0][j0][l1] + t1*prevDField[i0][j1][l1]) + s1*(t0*prevDField[i1][j0][l1] + t1*prevDField[i1][j1][l1]);

					dField[i][j][l] = u0*temp0 + u1*temp1;
					*/

					// This should be the proper computation

					float temp0 = s0 * (t0 * (u0 * prevDField[i0][j0][l0] +
											u1 * prevDField[i0][j0][l1])
										+(t1 * (u0 * prevDField[i0][j1][l0] +
											u1 * prevDField[i0][j1][l1])));

					float temp1 = s1 * (t0 * (u0 * prevDField[i1][j0][l0] +
											u1 * prevDField[i1][j0][l1])
										+(t1 * (u0 * prevDField[i1][j1][l0] +
											u1 * prevDField[i1][j1][l1])));

					dField[i][j][l] = temp0 + temp1;

				}
			}
		}
		setBoundary(N, b, dField);
	}

	void FluidSolver3D::project (int N, std::vector<std::vector<std::vector<float> > > & uField, std::vector<std::vector<std::vector<float> > > & vField, std::vector<std::vector<std::vector<float> > > & wField,  std::vector<std::vector<std::vector<float> > > & pField, std::vector<std::vector<std::vector<float> > > & divField) {

		float scale = -0.5f;

		for (unsigned int i = 1; i <= N; i++) { 
				for (unsigned int j = 1; j <= N; j++) {
					for (unsigned int l = 1; l <= N; l++) {
						divField[i][j][l] = scale * (uField[i+1][j][l] - uField[i-1][j][l] + vField[i][j+1][l] - vField[i][j-1][l] + wField[i][j][l+1] - wField[i][j][l-1]) / N;
						pField[i][j][l] = 0.0f;
					}
				}
		}
			
		setBoundary(N, 0, divField);
		setBoundary(N, 0, pField);

		// TO TAKE INTO ACCOUNT 6-NEIGHBOURHOOD!
		linearSolve(N, 0, pField, divField, 1, 6);

		scale = 0.5f;

		for (unsigned int i = 1; i <= N; i++) { 
			for (unsigned int j = 1; j <= N; j++) {
				for (unsigned int l = 1; l <= N; l++) {
					uField[i][j][l] -= scale * (pField[i+1][j][l] - pField[i-1][j][l]) * N;
					vField[i][j][l] -= scale * (pField[i][j+1][l] - pField[i][j-1][l]) * N;
					wField[i][j][l] -= scale * (pField[i][j][l+1] - pField[i][j][l-1]) * N;
				}
			}
		}

		setBoundary(N, 1, uField);
		setBoundary(N, 2, vField);
		setBoundary(N, 3, wField);
	}

	/* The vorticity confinement for each cell is computed to preserve turbulences that
	 * are suppressed due to artificial numerical dissipation.
	 * The forces at points (i, j, l) are given by N x omega, where N are normalized
	 * vorticity location vectors pointing from lower vorticity concentrations to
	 * higher vorticity concentrations and omega is the vorticity.
	 * The forces are the scaled by multiplying by epsilon, which gives control
	 * on amount of small scale detail to be added back in the flow
	 */
	void FluidSolver3D::vorticityConfinement (float dt) {

		// Uses the prev velocity fields as temporary buffers, since
		// the computation happens before the swap
		
		float dt0 = dt * vorticityConfinementEps;

		float x, y, z;
		
		for (unsigned int i = 1; i < N; i++) { 
			for (unsigned int j = 1; j < N; j++) {
				for (unsigned int l = 1; l < N; l++) {

					x = _prevUField[i][j][l] = (_wField[i][j+1][l] - _wField[i][j-1][l]) * 0.5f - (_vField[i][j][l+1] - _vField[i][j][l-1]) * 0.5f;

					y = _prevVField[i][j][l] = (_uField[i][j][l+1] - _uField[i][j][l-1]) * 0.5f - (_wField[i+1][j][l] - _wField[i-1][j][l]) * 0.5f;

					z = _prevWField[i][j][l] = (_vField[i+1][j][l] - _vField[i-1][j][l]) * 0.5f - (_uField[i][j+1][l] - _uField[i][j-1][l]) * 0.5f;

					_vorticityField[i][j][l] = sqrtf(x*x+y*y+z*z);
				}
			}
		}

		for (unsigned int i = 1; i < N; i++) { 
			for (unsigned int j = 1; j < N; j++) {
				for (unsigned int l = 1; l < N; l++) {
					float Nx = (_vorticityField[i+1][j][l] - _vorticityField[i-1][j][l]) * 0.5f; 
					float Ny = (_vorticityField[i][j+1][l] - _vorticityField[i][j-1][l]) * 0.5f;
					float Nz = (_vorticityField[i][j][l+1] - _vorticityField[i][j][l-1]) * 0.5f;
					// add small epsilon to avoid / 0
					float lenght = 1.0f/(sqrtf(Nx*Nx+Ny*Ny+Nz*Nz)+0.0000001f);
					
					Nx *= lenght;
					Ny *= lenght;
					Nz *= lenght;
					

					/* Uses values from the temporary fields */
					_uField[i][j][l] += (Ny * _prevWField[i][j][l] - Nz * _prevVField[i][j][l]) * dt0;
					_vField[i][j][l] += (Nz * _prevUField[i][j][l] - Nx * _prevWField[i][j][l]) * dt0;
					_wField[i][j][l] += (Nx * _prevVField[i][j][l] - Ny * _prevUField[i][j][l]) * dt0;



				}
			}
		}

	}

	/* Implements the buoyancy force computation, as described in "Visual simulation
	 * of smoke" by Stam et al., with two additional assumptions.
	 * 1) Temperature is assumed being a fraction of density at each point.
	 * (This is in general not true. In fact, assuming an ideal gas at constant pressure
	 * the temperature would be inversely proportional to density). Hence the
	 * implementation is simplified a bit, with no need to keep and
	 * update an additional field to track temperature at each voxel.
	 * 2) Ambient temperature is assumed constant over all grid.
	 * The formula to compute the force is:
	 * -alpha rho[i,j,l] + beta (temp[i,j,l] - ambientTemperature)
	 * The result is then applied to the v velocity field simply multiplying the
	 * force by dt timestep parameter, as described in the paper.
	 */
	void FluidSolver3D::addBuoyancy (std::vector<std::vector<std::vector<float> > >& densityField, std::vector<std::vector<std::vector<float> > >& vField, float dt) {
		
		for (unsigned int i = 0; i < size; i++)
			for (unsigned int j = 0; j < size; j++)
				for (unsigned int l = 0; l < size; l++) {
					float temp = 0.6 * densityField[i][j][l];
					vField[i][j][l] += ((-buoyancyAlphaFactor * densityField[i][j][l]) + buoyancyBetaFactor * (temp - ambientTemp)) * dt;
				}

	}

	void FluidSolver3D::densityStep () {
		addSource(_prevDensityField, dt, _densityField);
		_prevDensityField.swap(_densityField);
		diffuse(N, 0, _densityField, _prevDensityField, diff, dt);
		_prevDensityField.swap(_densityField);
		advect(N, 0, _densityField, _prevDensityField, _uField, _vField, _wField, dt);
	}

	void FluidSolver3D::velocityStep () {
		
		addSource(_prevUField, dt, _uField); 
		addSource(_prevVField, dt, _vField);
		addSource(_prevWField, dt, _wField);

		/**********************************
		* VORTICITY CONFINEMENT EXTENSION *
		**********************************/
		if (confineVorticity) {
			addBuoyancy(_densityField, _vField, dt);
			vorticityConfinement(dt);
		}
	
		_prevUField.swap(_uField);
		diffuse(N, 1, _uField, _prevUField, visc, dt);

		_prevVField.swap(_vField);
		diffuse(N, 2, _vField, _prevVField, visc, dt);

		_prevWField.swap(_wField);
		diffuse(N, 3, _wField, _prevWField, visc, dt);

		project(N, _uField, _vField, _wField, _prevUField, _prevVField);

		_prevUField.swap(_uField);
		_prevVField.swap(_vField);
		_prevWField.swap(_wField);
		
		advect(N, 1, _uField, _prevUField, _prevUField, _prevVField, _prevWField, dt);
		advect(N, 2, _vField, _prevVField, _prevUField, _prevVField, _prevWField, dt);
		advect(N, 3, _wField, _prevWField, _prevUField, _prevVField, _prevWField, dt);
		project(N, _uField, _vField, _wField, _prevUField, _prevVField);
	}

	// create a volume texture with n^3 texels and base radius r
	GLuint FluidSolver3D::get3DTexture () {
		int n = N + 2;
		float r = 0.5f; 

		GLuint texid;
		glGenTextures(1, &texid);

		GLenum target = GL_TEXTURE_3D;
		GLenum filter = GL_LINEAR;
		GLenum address = GL_CLAMP_TO_BORDER;

		glBindTexture(target, texid);

		glTexParameteri(target, GL_TEXTURE_MAG_FILTER, filter);
		glTexParameteri(target, GL_TEXTURE_MIN_FILTER, filter);

		glTexParameteri(target, GL_TEXTURE_WRAP_S, address);
		glTexParameteri(target, GL_TEXTURE_WRAP_T, address);
		glTexParameteri(target, GL_TEXTURE_WRAP_R, address);

		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

		GLbyte* data = new GLbyte[n*n*n];
		GLbyte* ptr = data;

		//float frequency = 3.0f / n;
		//float center = n / 2.0f + 0.5f;

		//for(int x=0; x < n; x++) {
		//	for (int y=0; y < n; ++y) {
		//		for (int z=0; z < n; ++z) {
		//			float dx = center-x;
		//			float dy = center-y;
		//			float dz = center-z;

		//			float randomOffset = RTA::RandomFloat(0.0f, 1.0f);
		//			//float off = fabsf(Perlin3D(x*frequency,
		//			//	y*frequency,
		//			//	z*frequency,
		//			//	5,
		//			//	0.5f));

		//			float d = sqrtf(dx*dx+dy*dy+dz*dz)/(n);

		//			*ptr++ = ((d-randomOffset) < r)?255:0;
		//		}
		//	}
		//}

		for (unsigned int l = 0; l < _densityField.size(); l++)
			for (unsigned int j = 0; j < _densityField.size(); j++)
				for (unsigned int i = 0; i < _densityField.size(); i++) {
					// forces a discrete value for the density
					// Could it be done better by using 4 bytes for each float
					//*ptr++ = (GLbyte)((int)((_densityField[i][j][l]) * (256.0 / 3.0)));
						
					
					
					//std::cout << "Value int in 0 255: " << intConv << std::endl;
					//	GLbyte byteConv = (GLbyte)intConv;
					//if (_densityField[i][j][l] > 1.0f) {
					//	std::cout << "Value: " << _densityField[i][j][l] << std::endl;
					//	int intConv = (int)((_densityField[i][j][l]) * (256.0 / 2.0));
					//	std::cout << "Value int in 0 255: " << intConv << std::endl;
					//	GLbyte byteConv = (GLbyte)intConv;
					//	std::cout << "Value byte: " << byteConv << std::endl;
					//}
					//*ptr++ = _densityField[i][j][l];
				}
		// upload
		glTexImage3D(target,
			0,
			GL_LUMINANCE,
			n,
			n,
			n,
			0,
			GL_LUMINANCE,
			GL_UNSIGNED_BYTE,
			data);

		delete[] data;

		return texid;
	}

	void FluidSolver3D::drawContainer () {

		float x, y, z;

		float dimOffset = N/2.0f;
		x = origin.x - dimOffset;
		y = origin.y - dimOffset;
		z = origin.z - dimOffset;

		float h = N + 2;

		glColor3f(0.0f, 1.0f, 1.0f);

		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glLineWidth(1.5f);
		glBegin(GL_QUADS);


		// bottom
		glVertex3f(x, y, z);
		glVertex3f(x+h, y, z);
		glVertex3f(x+h, y, z+h);
		glVertex3f(x, y, z+h);

		// front
		glVertex3f(x, y, z);
		glVertex3f(x+h, y, z);
		glVertex3f(x+h, y+h, z);
		glVertex3f(x, y+h, z);

		// right
		glVertex3f(x+h, y, z);
		glVertex3f(x+h, y, z+h);
		glVertex3f(x+h, y+h, z+h);
		glVertex3f(x+h, y+h, z);

		// left
		glVertex3f(x, y, z);
		glVertex3f(x, y, z+h);
		glVertex3f(x, y+h, z+h);
		glVertex3f(x, y+h, z);

		// back
		glVertex3f(x, y, z+h);
		glVertex3f(x+h, y, z+h);
		glVertex3f(x+h, y+h, z+h);
		glVertex3f(x, y+h, z+h);

		// top
		glVertex3f(x, y+h, z);
		glVertex3f(x+h, y+h, z);
		glVertex3f(x+h, y+h, z+h);
		glVertex3f(x, y+h, z+h);

		glEnd ();
		glLineWidth(1.0f);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		glColor3f(1.0f, 1.0f, 1.0f);
	}

	// TODO implement with 3D textures
	void FluidSolver3D::drawDensity () {
		static bool init = false;
		if (!init) {
			textureID = get3DTexture();
			std::cout << "3D texture initialized" << std::endl;
			init = true;
		}

		glEnable(GL_TEXTURE_3D);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_3D, textureID);

		float x, y, z;

		float dimOffset = N/2.0f;
		RTA::Vector3 offset = RTA::Vector3(origin.x - dimOffset, origin.y - dimOffset, origin.z - dimOffset);
		x = offset.x;
		y = offset.y;
		z = offset.z;
		float h = (float)(N + 2);

		glColor3f(1.0f,1.0f,1.0f); 

		glBegin(GL_QUADS);
		// bottom
		glTexCoord3f(0.0f, 0.0f, 0.0f); glVertex3f(x, y, z);

		glTexCoord3f(1.0f, 0.0f, 0.0f); glVertex3f(x+h, y, z);

		glTexCoord3f(1.0f, 0.0f, 1.0f); glVertex3f(x+h, y, z+h);

		glTexCoord3f(0.0f, 0.0f, 1.0f); glVertex3f(x, y, z+h);

		// front
		glTexCoord3f(0.0f, 0.0f, 0.0f); glVertex3f(x, y, z);

		glTexCoord3f(1.0f, 0.0f, 0.0f); glVertex3f(x+h, y, z);

		glTexCoord3f(1.0f, 1.0f, 0.0f); glVertex3f(x+h, y+h, z);

		glTexCoord3f(0.0f, 1.0f, 0.0f); glVertex3f(x, y+h, z);

		// right
		glTexCoord3f(1.0f, 0.0f, 0.0f); glVertex3f(x+h, y, z);

		glTexCoord3f(1.0f, 0.0f, 1.0f); glVertex3f(x+h, y, z+h);

		glTexCoord3f(1.0f, 1.0f, 1.0f); glVertex3f(x+h, y+h, z+h);

		glTexCoord3f(1.0f, 1.0f, 0.0f); glVertex3f(x+h, y+h, z);

		// left
		glTexCoord3f(0.0f, 0.0f, 0.0f); glVertex3f(x, y, z);

		glTexCoord3f(0.0f, 0.0f, 1.0f); glVertex3f(x, y, z+h);

		glTexCoord3f(0.0f, 1.0f, 1.0f); glVertex3f(x, y+h, z+h);

		glTexCoord3f(0.0f, 1.0f, 0.0f); glVertex3f(x, y+h, z);

		// back
		glTexCoord3f(0.0f, 0.0f, 1.0f); glVertex3f(x, y, z+h);

		glTexCoord3f(1.0f, 0.0f, 1.0f); glVertex3f(x+h, y, z+h);

		glTexCoord3f(1.0f, 1.0f, 1.0f); glVertex3f(x+h, y+h, z+h);

		glTexCoord3f(0.0f, 1.0f, 1.0f); glVertex3f(x, y+h, z+h);

		// top
		glTexCoord3f(0.0f, 1.0f, 0.0f); glVertex3f(x, y+h, z);

		glTexCoord3f(1.0f, 1.0f, 0.0f); glVertex3f(x+h, y+h, z);

		glTexCoord3f(1.0f, 1.0f, 1.0f); glVertex3f(x+h, y+h, z+h);

		glTexCoord3f(0.0f, 1.0f, 1.0f); glVertex3f(x, y+h, z+h);

		glEnd ();

		glDisable(GL_TEXTURE_3D);
		glColor3f(1.0f, 1.0f, 1.0f);

	}

	void FluidSolver3D::drawDensityWithBoxes () {
		float x, y, z;
		float d000, d001, d010, d100, d011, d101, d110, d111;

		float dimOffset = N/2.0f;
		RTA::Vector3 offset = RTA::Vector3(origin.x - dimOffset, origin.y - dimOffset, origin.z - dimOffset);

		float h = 1.0;

		// Using different increases in brightness produces better volumetric rendering
		float brightnessBase1 = 0.0f;
		float brightnessBase2 = 0.0f;	// allows better volume rendering by applying different color bases to the vertices

		glBegin(GL_QUADS);

			for (unsigned int i = 0; i <= N; i++) {
				x = offset.x + i;
				for (unsigned int j = 0; j <= N; j++) {
					y = offset.y + j;
					for (unsigned int l = 0; l <= N; l++) {
						//if (_densityField[i][j][l] > 0.0f) {
							z = offset.z + l;

							float alpha = _densityField[i][j][l];
						/*	d000 = brightnessBase1 + _densityField[i][j][l];
							d001 = brightnessBase2 + _densityField[i][j][l];
							d010 = brightnessBase1 + _densityField[i][j][l];
							d100 = brightnessBase2 + _densityField[i][j][l];
							d011 = brightnessBase1 + _densityField[i][j][l];
							d101 = brightnessBase2 + _densityField[i][j][l];
							d110 = brightnessBase1 + _densityField[i][j][l];
							d111 = brightnessBase2 + _densityField[i][j][l];*/
							d000 = brightnessBase1 + _densityField[i][j][l];
							d001 = brightnessBase2 + _densityField[i][j][l+1];
							d010 = brightnessBase1 + _densityField[i][j+1][l];
							d100 = brightnessBase2 + _densityField[i+1][j][l];
							d011 = brightnessBase1 + _densityField[i][j+1][l+1];
							d101 = brightnessBase2 + _densityField[i+1][j][l+1];
							d110 = brightnessBase1 + _densityField[i+1][j+1][l];
							d111 = brightnessBase2 + _densityField[i+1][j+1][l+1];

							////////////////////////////////////////
							/// HACK: Drawing only the front face //
							////////////////////////////////////////
							/* Proper drawing with alpha blending would require
							 * to (re)sort the density cubes based on camera position
							 * and orientation, and drawing back to front the
							 * faces of the cubes visible from that point
							 */

							// bottom
							/*glColor4f(d000, d000, d000, alpha); 
							glVertex3f(x, y, z);

							glColor4f(d100, d100, d100, alpha);
							glVertex3f(x+h, y, z);

							glColor4f(d101, d101, d101, alpha);
							glVertex3f(x+h, y, z+h);

							glColor4f(d001, d001, d001, alpha);
							glVertex3f(x, y, z+h);*/

							// front
							glColor4f(d000, d000, d000, alpha); 
							glVertex3f(x, y, z);

							glColor4f(d100, d100, d100, alpha);
							glVertex3f(x+h, y, z);

							glColor4f(d110, d110, d110, alpha);
							glVertex3f(x+h, y+h, z);

							glColor4f(d010, d010, d010, alpha);
							glVertex3f(x, y+h, z);

							//// right
							//glColor4f(d100, d100, d100, alpha); 
							//glVertex3f(x+h, y, z);

							//glColor4f(d101, d101, d101, alpha);
							//glVertex3f(x+h, y, z+h);

							//glColor4f(d111, d111, d111, alpha);
							//glVertex3f(x+h, y+h, z+h);

							//glColor4f(d110, d110, d110, alpha);
							//glVertex3f(x+h, y+h, z);

							//// left
							//glColor4f(d000, d000, d000, alpha); 
							//glVertex3f(x, y, z);

							//glColor4f(d001, d001, d001, alpha);
							//glVertex3f(x, y, z+h);

							//glColor4f(d011, d011, d011, alpha);
							//glVertex3f(x, y+h, z+h);

							//glColor4f(d010, d010, d010, alpha);
							//glVertex3f(x, y+h, z);

							//// back
							//glColor4f(d001, d001, d001, alpha);
							//glVertex3f(x, y, z+h);

							//glColor4f(d101, d101, d101, alpha);
							//glVertex3f(x+h, y, z+h);

							//glColor4f(d111, d111, d111, alpha);
							//glVertex3f(x+h, y+h, z+h);

							//glColor4f(d011, d011, d011, alpha);
							//glVertex3f(x, y+h, z+h);

							//// top
							//glColor4f(d010, d010, d010, alpha);
							//glVertex3f(x, y+h, z);

							//glColor4f(d110, d110, d110, alpha);
							//glVertex3f(x+h, y+h, z);

							//glColor4f(d111, d111, d111, alpha);
							//glVertex3f(x+h, y+h, z+h);

							//glColor4f(d011, d011, d011, alpha);
							//glVertex3f(x, y+h, z+h);




							//float red = l == N ? 1.0f : 0.0f;
							//// bottom
							//glColor4f(red, 0.0f, 0.0f, alpha); 
							//glVertex3f(x, y, z);

							//glColor4f(red, 0.0f, 0.0f, alpha); 
							//glVertex3f(x+h, y, z);

							//glColor4f(red, 0.0f, 0.0f, alpha); 
							//glVertex3f(x+h, y, z+h);

							//glColor4f(red, 0.0f, 0.0f, alpha); 
							//glVertex3f(x, y, z+h);

							//// front
							//glColor4f(red, 0.0f, 0.0f, alpha); 
							//glVertex3f(x, y, z);

							//glColor4f(red, 0.0f, 0.0f, alpha); 
							//glVertex3f(x+h, y, z);

							//glColor4f(red, 0.0f, 0.0f, alpha); 
							//glVertex3f(x+h, y+h, z);

							//glColor4f(red, 0.0f, 0.0f, alpha); 
							//glVertex3f(x, y+h, z);

							//// right
							//glColor4f(red, 0.0f, 0.0f, alpha); 
							//glVertex3f(x+h, y, z);

							//glColor4f(red, 0.0f, 0.0f, alpha); 
							//glVertex3f(x+h, y, z+h);

							//glColor4f(red, 0.0f, 0.0f, alpha); 
							//glVertex3f(x+h, y+h, z+h);

							//glColor4f(red, 0.0f, 0.0f, alpha); 
							//glVertex3f(x+h, y+h, z);

							//// left
							//glColor4f(red, 0.0f, 0.0f, alpha); 
							//glVertex3f(x, y, z);

							//glColor4f(red, 0.0f, 0.0f, alpha); 
							//glVertex3f(x, y, z+h);

							//glColor4f(red, 0.0f, 0.0f, alpha); 
							//glVertex3f(x, y+h, z+h);

							//glColor4f(red, 0.0f, 0.0f, alpha); 
							//glVertex3f(x, y+h, z);

							//// back
							//glColor4f(red, 0.0f, 0.0f, alpha); 
							//glVertex3f(x, y, z+h);

							//glColor4f(red, 0.0f, 0.0f, alpha); 
							//glVertex3f(x+h, y, z+h);

							//glColor4f(red, 0.0f, 0.0f, alpha); 
							//glVertex3f(x+h, y+h, z+h);

							//glColor4f(red, 0.0f, 0.0f, alpha); 
							//glVertex3f(x, y+h, z+h);

							//// top
							//glColor4f(red, 0.0f, 0.0f, alpha); 
							//glVertex3f(x, y+h, z);

							//glColor4f(red, 0.0f, 0.0f, alpha); 
							//glVertex3f(x+h, y+h, z);

							//glColor4f(red, 0.0f, 0.0f, alpha); 
							//glVertex3f(x+h, y+h, z+h);

							//glColor4f(red, 0.0f, 0.0f, alpha); 
							//glVertex3f(x, y+h, z+h);
						//}
					}
				}
			}

		glEnd ();

		glColor4f(1.0f, 1.0f, 1.0f, 1.0f);

	}

	void FluidSolver3D::drawMarkerParticleSet () {
		
		float dimOffset = N/2.0f;
		RTA::Vector3 offset = RTA::Vector3(origin.x - dimOffset, origin.y - dimOffset, origin.z - dimOffset);

		glColor3f(1.0f, 0.0f, 0.0f);
		glPointSize(1.5f);

		glBegin (GL_POINTS);

			for (unsigned int i = 0; i < NUMOF3DMARKERPARTICLES; i++) {
				glVertex3f(offset.x + particles[i].x, offset.y + particles[i].y, offset.z + particles[i].z);
			}

		glEnd ();

		glColor3f(1.0f, 1.0f, 1.0f);
		glPointSize(1.0f);
	}


	void FluidSolver3D::drawVelocity () {
		float x, y, z;

		float dimOffset = N/2.0f;
		RTA::Vector3 offset = RTA::Vector3(origin.x - dimOffset, origin.y - dimOffset, origin.z - dimOffset);

		glLineWidth(1.0);
		
		glBegin (GL_LINES);

			glColor3f(1.0f, 1.0f, 0.0f);
		
			for (unsigned int i = 0; i <= N; i++) {
				x = offset.x + i;
				for (unsigned int j = 0; j <= N; j++) {
					y = offset.y + j;
					for (unsigned int l = 0; l <= N; l++) {
						z = offset.z + l;
						glVertex3f(x, y, z);
						glVertex3f(x+_uField[i][j][l], y+_vField[i][j][l], z + _wField[i][j][l]);
					}
				}
			}

			glColor3f(1.0f, 1.0f, 1.0f);

		glEnd ();

	}

}

