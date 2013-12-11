package net.maizegenetics.pal.phenotype;

public class NumericAttribute implements PhenotypeAttribute {
	private String name;
	private float[] values;
	private boolean[] missing;

	public float getFloatValue(int obs) {
		return values[obs];
	}
	
	public float[] getFloatValues() {
		return values;
	}
	
	@Override
	public Object getValue(int obs) {
		return new Float(values[obs]);
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
