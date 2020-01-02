package main;

public class CFDMain {
	
	private int cellNumberX=20;
	private int cellNumberY=20;
	
	private  double domainHeight=1;
	private  double domainLength=1;
	
	private double[] cellLeftPosX;
	private double[] cellBottomPosY;
	
	private double[] cellMidPosX;
	private double[] cellMidPosY;
	
	
	private int imin, imax;
	private int jmin, jmax;
	
	private double dx=0;
	private double dy=0;
	
	private double dxi=0;
	private double dyi=0;
	
	private double[][] uStar;
	private double[][] vStar;
	private double[][] u;
	private double[][] v;
	
	private double[][] p;
	private double[] pressureVector;
	
	private double density=1;
	
	private double dynamicVisc=1;
	private double kinematicVisc=dynamicVisc/density;
	
	private double[][] laplacianOperator;
	private double[] RightHandSide;
	private double[] pressureTerm;
	
	
	public CFDMain() {
		
		// Inititaliase cell position:
		cellLeftPosX 	= new double[cellNumberX+1];
		cellBottomPosY 	= new double[cellNumberY+1];
		cellMidPosX 	    = new double[cellNumberX+1];
		cellMidPosY     = new double[cellNumberY+1];
		
		cellLeftPosX[0] = 0;
		cellMidPosX[0]  = (domainLength/cellNumberX)/2;
		
		imin = 1 ;
		jmin = 1 ;
		imax = cellNumberX -1;
		jmax = cellNumberX -1;
		
		for(int i=1;i<cellNumberX-1;i++) {
			cellLeftPosX[i] = cellLeftPosX[i-1] + domainLength/(cellNumberX-1);
			cellMidPosX[i] = cellMidPosX[i-1] + (domainLength/cellNumberX);
			//System.out.println(cellMidPosX[i]);
		}
		cellBottomPosY[0] = 0;
		cellMidPosY[0]    = (domainHeight/cellNumberY)/2;
		for(int j=1;j<cellNumberY-1;j++) {
			cellBottomPosY[j] = cellBottomPosY[j-1] + domainHeight/(cellNumberY-1);
			cellMidPosY[j] = cellMidPosY[j-1] + (domainHeight/cellNumberY);
		}
		
		dx = cellLeftPosX[1] - cellLeftPosX[0];
		dy = cellBottomPosY[1] - cellBottomPosY[0];
		dxi=1/dx;
		dyi=1/dy;
		//System.out.println(""+dy);
		//---------------------------------------------------------------------------
		
		uStar = new double[cellNumberX+1][cellNumberY+1];
		vStar = new double[cellNumberX+1][cellNumberY+1];
		u	  = new double[cellNumberX+1][cellNumberY+1];
		v	  = new double[cellNumberX+1][cellNumberY+1];
		p	  = new double[cellNumberX+1][cellNumberY+1];
		pressureVector = new double[(cellNumberX+1)*(cellNumberY+1)];
		/**
		 * 
		 *  			Initialize uStar, Vstar, u and v 
		 * 
		 * 
		 */
		for(int i=0;i<cellNumberX-1;i++) {
			for(int j=0;j<cellNumberY-1;j++) {
				uStar[i][j] = 0;
				vStar[i][j] = 0;
				u[i][j] 		= 0;
				v[i][j] 		= 0;
				p[i][j] 		= 0;
			}
		}
		
		/**
		 * 
		 * 
		 * 			Initialize Laplacian Operator
		 * 
		 */
		int laplaceSize = cellNumberX * cellNumberY;
		laplacianOperator = new double[laplaceSize][laplaceSize];
		for(int i=0;i<laplaceSize;i++) {
			for(int j=0;j<laplaceSize;j++) {
				laplacianOperator[i][j]=0;
			}
		}
		
		
		/**
		 * 
		 *  			Set boundary condition
		 * 
		 * 
		 */
		double uu = 0;
		double ud = 0;
		double vl = 10;
		double vr = 0; 
		setBoundaryCondition( u,  v,  uu,  ud,  vl,  vr);
	}
	
	private void setBoundaryCondition(double[][] u, double[][] v, double uu, double ud, double vl, double vr) {
		for(int i=imin;i<imax;i++) {
			u[i][jmin-1] = u[i][jmin] - 2 * (u[i][jmin] - ud);
			u[i][jmax+1] = u[i][jmax] - 2 * (u[i][jmax] - uu);
		}
		for(int j=jmin;j<jmax;j++) {
			v[imin-1][j] = v[imin][j] - 2 * ( v[imin][j] - vl);
			v[imax+1][j] = v[imax][j] - 2 * ( v[imax][j] - vr);
		}
	}
	
	private void correctorStep(double[][] u, double[][] v, double[][] uStar, double[][] vStar, double dt) {
		for(int j=jmin;j<jmax;j++) {
			for(int i=imin;i<imax;i++) {
				u[i][j] = uStar[i][j] - dt/density * ( p[i][j] - p[i-1][j] ) * dxi ;
				v[i][j] = vStar[i][j] - dt/density * ( p[i][j] - p[i][j-1] ) * dyi ;
			}
		}
	}
	
	private double[][] computeUStar(double[][] u, double[][] v, double dt){
		double[][] uStar = new double[cellNumberX+1][cellNumberY+1];
		for(int i=1;i<cellNumberX;i++) {
			for(int j=1;j<cellNumberY-1;j++) {
				double vInterim = 0.25 * (v[i][j-1] +v[i-1][j+1] + v[i][j] + v[i][j+1] ) ;
				uStar[i][j] = u[i][j] + dt*( kinematicVisc * 
						( ( u[i-1][j] - 2 * u[i][j] + u[i+1][j] ) * dxi * dxi 
						+ ( u[i][j-1] - 2 * u[i][j] + u[i][j+1] ) * dyi * dyi )
						- u[i][j] * ( u[i+1][j] - u[i-1][j] ) * 0.5 * dxi 
						- vInterim * ( u[i][j+1] -u[i][j-1] ) * 0.5 * dyi );
				//System.out.println(uStar[i][j]);
			}
		}
		return uStar;
	}
	
	private double[][] computeVStar(double[][] u, double[][] v, double dt){
		double[][] vStar = new double[cellNumberX+1][cellNumberY+1];
		for(int i=1;i<cellNumberX;i++) {
			for(int j=1;j<cellNumberY-1;j++) {
				double uInterim = 0.25 * (u[i][j-1] +u[i+1][j-1] + u[i][j] + u[i+1][j] ) ;
				uStar[i][j] = v[i][j] + dt*( kinematicVisc * 
						( ( v[i-1][j] - 2 * v[i][j] + v[i+1][j] ) * dxi * dxi 
						+ ( v[i][j-1] - 2 * v[i][j] + v[i][j+1] ) * dyi * dyi )
						- uInterim * ( v[i+1][j] - v[i-1][j] ) * 0.5 * dxi 
						- v[i][j]  * ( v[i][j+1] - v[i][j-1] ) * 0.5 * dyi );
				//System.out.println(uStar[i][j]);
			}
		}
		return vStar;
	}
	
	private double[][] createLaplacianOperator(double[][] laplacianOperator){
		for(int j=0;j<cellNumberY;j++) {
			for(int i=0;i<cellNumberX;i++) {
				laplacianOperator[i+(j-1)*cellNumberX][i+(j-1)*cellNumberY] = 2 * dxi*dxi + 2 * dyi * dyi;
				for(int ii=i-1;ii<i+1;i=i+2) {
					if(ii>0 && ii<=cellNumberX) { // Interior point
						laplacianOperator[i+(j-1)*cellNumberX][ii+(j-1)*cellNumberX]=-dxi*dxi;
					} else {	 // Neuman conditions on boundary 
						laplacianOperator[i+(j-1)*cellNumberX][i+(j-1)*cellNumberX] = laplacianOperator[i+(j-1)*cellNumberX][i+(j-1)*cellNumberX]-dxi*dxi;
					}
				}
				for(int jj=j-1;jj<j+1;j=j+2) {
					if(jj>0 && jj<=cellNumberY) {
						laplacianOperator[i+(j-1)*cellNumberX][i+(jj-1)*cellNumberX] = -dyi * dyi;
					} else {
						laplacianOperator[i+(j-1)*cellNumberX][i+(j-1)*cellNumberX]  = laplacianOperator[i+(j-1)*cellNumberX][i+(j-1)*cellNumberX] = - dyi * dyi;
					}
				}
			}
		}
		laplacianOperator[0][0] = 1;
		return laplacianOperator;
	}
	
	private double[] createRightHandSide(double[] RightHandSide, double[][] uStar, double[][] vStar, double dt) {
		int n=0;
		for(int j=jmin;j<jmax;j++) {
			for(int i=imin;i<imax;i++) {
				RightHandSide[n] = -density/dt * ( (uStar[i+1][j] - uStar[i][j] ) * dxi + (vStar[i][j+1] - vStar[i][j] ) * dyi);
				n++;
			}
		}
		
		return RightHandSide;
	}
	
	private double[][] createPressure(double[] pressureTerm, double[][] p) {
		int n=0;
		for(int j=jmin; j<jmax;j++) {
			for(int i=imin;i<imax;i++) {
				p[i][j] = pressureTerm[n];
				n++;
			}
		}
		return p;
	}
	
	private void run() {
		double dt = 0.1;
		double tmax = 10;
		double t=0;
		
	
		
		laplacianOperator = createLaplacianOperator(laplacianOperator);
		p = createPressure(pressureTerm, p);
		
		while(t<tmax) {
			/**
			 * 	Predictor step
			 */
			 uStar = computeUStar(u, v, dt);
			 vStar = computeVStar(u, v, dt);
			/**
			 *  Forming right hand side of poisson equation 
			 */
			 RightHandSide = createRightHandSide(RightHandSide,  uStar, vStar,  dt);
			 
			 laplacianOperator = createLaplacianOperator(laplacianOperator);
			 /**
			  * 	Solving for pressure 
			  */
			 //pressureVector = L/R;
			 
			 p = createPressure( pressureTerm,  p);
			 /**
			  *  Corrector step 
			  */
			 correctorStep( u, v,  uStar,  vStar, dt);
		}
	}
	
	
	
	public static void main(String[] args) {
		
		CFDMain fluoFlow = new CFDMain();
		fluoFlow.run();
		
	}

}
