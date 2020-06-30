package mixedmodel;

import flanagan.roots.RealRootFunction;

public class RootSearchEmmadeltaREMLdLLwoZ implements RealRootFunction {
	private double[] lambda;
	private double[] etas;
	
//	emma.delta.REML.dLL.wo.Z <- function(logdelta, lambda, etas) {
//		  nq <- length(etas)
//		  delta <- exp(logdelta)
//		  etasq <- etas*etas
//		  ldelta <- lambda+delta
//		  return( 0.5*(nq*sum(etasq/(ldelta*ldelta))/sum(etasq/ldelta)-sum(1/ldelta)) )
//		}
	
	public double function(double logdelta) {
		double nq = this.etas.length;
		double delta = Math.exp(logdelta);
		double sum1=0,sum2=0,sum3=0;
		for(int i=0;i<etas.length;i++){
			  sum1+=etas[i]*etas[i]/((lambda[i]+delta)*(lambda[i]+delta));
			  sum2+=etas[i]*etas[i]/(lambda[i]+delta);
			  sum3+=1.0/(lambda[i]+delta);
		  }
		  return 0.5*(nq*sum1/sum2-sum3);
	}
	
	public void set_delta_REML_dLL_wo_Z( double[] lambda, double[] etas){
		this.lambda=lambda;
		this.etas=etas;
	}

}
