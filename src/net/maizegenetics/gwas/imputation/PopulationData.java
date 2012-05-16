package net.maizegenetics.gwas.imputation;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Pattern;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.util.BitSet;

public class PopulationData {
	public String name;
	public ArrayList<String> members;
	public BitSet snpIndex;
	public String parent1;
	public String parent2;
	public double contribution1;
	public double contribution2;
	public int Fgen;
	public double inbredCoef;
	public Alignment align;
	public byte[] alleleA;
	public byte[] alleleC;

	/**
	 * @param popFilename	the name of the file containing pedigree information for a group of populations
	 * @return	a HashMap of family names (keys) and associated PopulationData objects (values) containing information for the families in the pedigree file
	 */
	public static HashMap<String, PopulationData> readPedigreeFile(String popFilename) {
		Pattern tab = Pattern.compile("\t");
		HashMap<String, PopulationData> familyMap = new HashMap<String, PopulationData>(); 
		try {
			BufferedReader br = new BufferedReader(new FileReader(popFilename));
			br.readLine();
			String input;
			while ((input = br.readLine()) != null) {
				String[] info = tab.split(input);
				PopulationData family = familyMap.get(info[0]);
				if (family == null) {
					family = new PopulationData ();
					family.members = new ArrayList<String>();
					family.members.add(info[2]);  //add parents to family members
					family.members.add(info[3]);
					family.members.add(info[1]);
					family.parent1 = info[2];
					family.parent2 = info[3];
					family.contribution1 = Double.parseDouble(info[4]);
					family.contribution2 = Double.parseDouble(info[5]);
					familyMap.put(info[0], family);
				}
				else family.members.add(info[1]);
			}
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return familyMap;
	}
}