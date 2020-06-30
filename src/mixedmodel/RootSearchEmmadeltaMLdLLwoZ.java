package mixedmodel;

import flanagan.roots.RealRootDerivFunction;
import flanagan.roots.RealRootFunction;

public class RootSearchEmmadeltaMLdLLwoZ implements RealRootFunction {
	
//	private double a = 1.0D; 				// sample code from flanagan
//    private double b = 2.0D; 				// sample code from flanagan
    
    private double[] lambda;
    private double[] etas;
    private double[] xi;
    
//    public double[ ] function(double x){  	// sample code from flanagan
//             double[ ] y = new double[2]; 	// sample code from flanagan
//
//             y[0] = a*Math.pow(x,b) - 2.0D;	// sample code from flanagan
//             y[1] = a*b*Math.pow(x,b-1.0D);	// sample code from flanagan
//             return y;						// sample code from flanagan
//    } 
//
//    public void setA(double a){				// sample code from flanagan
//             this.a = a;					// sample code from flanagan
//    } 
//
//    public void setB(double b){				// sample code from flanagan
//             this.b = b;					// sample code from flanagan
//    }
	
    /*
	 *  translated from EMMA/emma.delta.ML.dLL.wo.Z <- function(logdelta, lambda, etas, xi)
	 */
	public double function(double logdelta) {
		  double n = xi.length;
		  double delta = Math.exp(logdelta);
		  //etasq <- etas*etas
		  //ldelta <- lambda+delta
		  double sum1=0,sum2=0,sum3=0;
		  for(int i=0;i<etas.length;i++){
			  sum1+=etas[i]*etas[i]/((lambda[i]+delta)*(lambda[i]+delta));
			  sum2+=etas[i]*etas[i]/(lambda[i]+delta);
		  }for(int i=0;i<xi.length;i++){
			  sum3+=1.0/(xi[i]+delta);
		  }
		  return 0.5*(n*sum1/sum2-sum3);
		}
	
	public void set_delta_ML_dLL_wo_Z( double[] lambda, double[] etas, double[] xi){
		this.lambda=lambda;
		this.etas=etas;
		this.xi=xi;
	}

}
