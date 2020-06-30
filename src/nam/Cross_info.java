package nam;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashMap;
import java.util.Iterator;


public class Cross_info {
	public HashMap<String, String> ecoID_natureName;
	public HashMap<String, String> natureName_ecoID;
	public HashMap<String, String> ecoID_crossid;
	
	public Cross_info(String cross_id){
		this.ecoID_natureName = new HashMap<String,String>();
		this.natureName_ecoID = new HashMap<String, String>();
		this.ecoID_crossid = new HashMap<String, String>();
		try{
//			 read founders' ecotype ID, put into hashmap for ecotype ID and nature name
			 BufferedReader fid_br = new BufferedReader(new FileReader(cross_id));
			 String lfid = fid_br.readLine();
			 while(lfid!=null){
				 String[] array = lfid.split("\t");
				 this.ecoID_natureName.put(array[1], array[0]);
				 this.ecoID_natureName.put(array[3], array[2]);
				 this.natureName_ecoID.put(array[0], array[1]);
				 this.natureName_ecoID.put(array[2], array[3]);
				 this.ecoID_crossid.put(array[1]+"_"+array[3], array[4]);
				 this.ecoID_crossid.put(array[3]+"_"+array[1], array[4]);				 
				 lfid = fid_br.readLine();
			 }
		}catch(Exception e){e.printStackTrace();}	
	}

}
