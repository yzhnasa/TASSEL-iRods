// RateDistribution.java
//
// (c) 1999-2000 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package net.maizegenetics.pal.substmodel;

import net.maizegenetics.pal.io.FormattedOutput;
import net.maizegenetics.pal.report.Report;
import net.maizegenetics.pal.util.Parameterized;

import java.io.PrintWriter;
import java.io.Serializable;


/**
 * abstract base class for models of rate variation over sites
 * employing a discrete rate distribution
 *
 * @version $Id: RateDistribution.java,v 1.1 2007/01/12 03:26:16 tcasstevens Exp $
 *
 * @author Korbinian Strimmer
 * @author Alexei Drummond
 */

public abstract class RateDistribution //extends PalObjectListener.EventGenerator
        implements Parameterized, Report, Cloneable, Serializable
{
	//
	// Public stuff
	//

	/** number of rate categories*/
	public int numRates;

	/** rates of each rate category */
	public double[] rate;

	/** probability of each rate */
	public double[] probability;

	//
	// Protected stuff
	//

	protected FormattedOutput format;

	//
	// Serialization code
	//
	private static final long serialVersionUID= -5584969247361304141L;

	//serialver -classpath ./classes net.maizegenetics.pal.substmodel.RateDistribution
	private void writeObject(java.io.ObjectOutputStream out) throws java.io.IOException {
		out.writeByte(1); //Version number
		out.writeObject(rate);
		out.writeObject(probability);
	}

	private void readObject(java.io.ObjectInputStream in) throws java.io.IOException, ClassNotFoundException{
		byte version = in.readByte();
		switch(version) {
			default : {
				rate= (double[])in.readObject();
				probability= (double[])in.readObject();
				numRates= rate.length;
				format = FormattedOutput.getInstance();
				break;
			}
		}
	}

	public final int getNumberOfRates() {
		return numRates;
	}

	/**
	 * construct discrete distribution
	 *
	 *  @param n number of rate categories
	 */
	public RateDistribution(int n)
	{
		format = FormattedOutput.getInstance();

		numRates = n;
		rate = new double[n];
		probability = new double[n];
	}

	// interface Report (remains abstract)

	// interface Parameterized (remains abstract)

	protected void printRates(PrintWriter out)
	{
		out.println("Relative rates and their probabilities:\n");
		format.displayIntegerWhite(out, numRates);
		out.println("   Rate      Probability");
		for (int i = 0; i < numRates; i++)
		{
			format.displayInteger(out, i+1, numRates);
			out.print("   ");
			format.displayDecimal(out, rate[i], 5);
			out.print("   ");
			format.displayDecimal(out, probability[i], 5);
			out.println();
		}
	}
	/**
	 * The non direct access method
	 */
	public final double[] getCategoryProbabilities() { return probability;	}
	public Object clone() {
		try {
			RateDistribution rd = (RateDistribution)super.clone();
			return rd;
		} catch (CloneNotSupportedException e) {
			// this shouldn't happen
			throw new InternalError();
		}
	}

}
