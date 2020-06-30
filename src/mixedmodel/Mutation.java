package mixedmodel;

public class Mutation {
	
	public int chr;
	public int loc;
	public String ref;
	public String alt;
	public double score;
	
	public boolean compound_het;
	public int chr2;
	public int loc2;
	public String ref2;
	public String alt2;
	public double score2;
	
	public Mutation(int chr, int loc, String ref, String alt, double score){
		this.chr=chr;
		this.loc=loc;
		this.ref= ref;
		this.alt= alt;
		this.score= score;
		this.compound_het=false;
	}
	
	public Mutation(int chr, int loc, String ref, String alt, double score,
			int chr2, int loc2, String ref2, String alt2, double score2){
		this.chr=chr;
		this.loc=loc;
		this.ref= ref;
		this.alt= alt;
		this.score= score;
		
		this.compound_het=true;
		this.chr2=chr2;
		this.loc2=loc2;
		this.ref2= ref2;
		this.alt2= alt2;
		this.score2= score2;
	}
	
	public Mutation(Mutation m1, Mutation m2){
		if(m1.compound_het||m2.compound_het){
			System.out.println("Error: can't combine compound hets into another compound het!!!");
			return;
		};
		this.chr=m1.chr;
		this.loc=m1.loc;
		this.ref= m1.ref;
		this.alt= m1.alt;
		this.score= m1.score;
		
		this.compound_het=true;
		this.chr2=m2.chr;
		this.loc2=m2.loc;
		this.ref2= m2.ref;
		this.alt2= m2.alt;
		this.score2= m2.score;
	}
	
	public String output(){		
		String out=this.chr+":"+this.loc+":"+this.ref+":"+this.alt+":"+this.score;
		if(compound_het){
			out=out+"/"+this.chr2+":"+this.loc2+":"+this.ref2+":"+this.alt2
					+":"+this.score2;
		}return out;
	}
	
}
