package net.maizegenetics.pal.phenotype;

import com.google.common.collect.ImmutableBiMap;

public class CategoricalAttribute implements PhenotypeAttribute {

	private String name;
	private int[] values;
	private ImmutableBiMap<String, Integer> labelBimap;
	private boolean[] missing;

	public int getIntegerValue(int obs) {
		return values[obs];
	}
	
	public int[] getIntegerValues() {
		return values;
	}
	
	String getAttrLabel(int index) {
		return labelBimap.inverse().get(index);
	}
	
	int getIndexForAttrLabel(String label) {
		return labelBimap.get(label);
	}

	String[] getValuesAsLabels() {
		int n = values.length;
		String[] labels = new String[n];
		ImmutableBiMap<Integer, String> reverseMap = labelBimap.inverse();
		for (int i = 0; i < n; i++) {
			labels[i] = reverseMap.get(i);
		}
		return labels;
	}
	
	@Override
	public Object getValue(int obs) {
		return labelBimap.inverse().get(values[obs]);
	}

	@Override
	public Object getValues() {
		return values;
	}

	@Override
	public boolean isMissing(int obs) {
		return missing[obs];
	}

	@Override
	public boolean[] getMissing() {
		return missing;
	}

	@Override
	public String getName() {
		return name;
	}

	@Override
	public int getSize() {
		return values.length;
	}

}
