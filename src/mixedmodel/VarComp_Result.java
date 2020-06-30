package mixedmodel;

import java.io.BufferedWriter;

public class VarComp_Result {

    public double delta;
    public double sigma_g;
    public double sigma_e;
    public double ml;
    public double[] beta;
    public double pvalue;
    public double R2;

    public void write2file(BufferedWriter bw) {
        try {
            //pvalue,AdjustedR2,coefficient,sigma_g2
            bw.write(this.pvalue + "," + this.R2 + "," + this.beta[1] + "," + this.sigma_g);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public String toString() {
        return (this.pvalue + "," + this.R2 + "," + this.beta[1] + "," + this.sigma_g);
    }

}
