package net.maizegenetics.pal.util;
import java.io.Serializable;
/**
 * Created using IntelliJ IDEA.
 * Author: Peter Bradbury
 * Date: May 18, 2004
 * Time: 4:04:46 PM
 */
//todo should this be moved into net.maizegenetics.pal
public class DataDefinition implements Serializable {
    //fields
    private int NumberOfHeaders;
    private int NumberOfTraits;
    private String[] headerLabel;
    private boolean[] headerUse;
    private boolean[] headerIsTrait;
    private String[] traitType;

    public static final String DATA = "Data";
    public static final String FACTOR = "Factor";
    public static final String COVAR= "Covariate";
    public static final String EXCLUDE = "Exclude";

    //constructors
    public DataDefinition(int NumberOfHeaders, int NumberOfTraits) {
        setNumberOfHeaders(NumberOfHeaders);
        setNumberOfTraits(NumberOfTraits);
        if (NumberOfHeaders > 0 ) headerLabel[0] = "Trait";
        if (NumberOfHeaders > 1) headerLabel[1] = "env";
    }

    //methods
    public String getHeaderName(int headerNumber) {
        return headerLabel[headerNumber];
    }

    public String getUseableHeaderName(int headerNumber) {
        int counter = -1;
        for (int i = 0; i < NumberOfHeaders; i++) {
            if (headerUse[i]) counter++;
            if (counter == headerNumber) return headerLabel[i];
        }
        return "";
    }

    public int getNumberOfHeaders() {
        return NumberOfHeaders;
    }

    public int getNumberOfUseableHeaders() {
        int n = 0;
        for (int i = 0; i < NumberOfHeaders; i++) {
            if (headerUse[i]) n++;
        }

        return n;
    }
    public int getNumberOfTraits() {
        return NumberOfTraits;
    }

    public String getTraitType(int traitNumber) {
        return traitType[traitNumber];
    }

    public boolean isHeaderUseable(int headerNumber) {
        return headerUse[headerNumber];
    }

    public boolean isHeaderTrait(int headerNumber) {
        return headerIsTrait[headerNumber];
    }

    public void setHeaderLabel(int headerNumber, String name) {
        headerLabel[headerNumber] = name;
    }

    public void useHeader(int headerNumber, boolean use) {
        headerUse[headerNumber] = use;
    }

    public void setIsTrait(int headerNumber, boolean isTrait) {
        headerIsTrait[headerNumber] = isTrait;
    }

    public void setTraitType(int traitNumber, String type) {
        traitType[traitNumber] = type;
    }

    protected void setNumberOfHeaders(int nHeaders) {
        this.NumberOfHeaders = nHeaders;
        headerLabel = new String[NumberOfHeaders];
        headerUse = new boolean[NumberOfHeaders];
        headerIsTrait = new boolean[NumberOfHeaders];
        setIsTrait(0, true);
        useHeader(0,true);
        if (NumberOfHeaders > 1) {
            headerLabel[1] = "env";
        }
        for (int i = 1; i < NumberOfHeaders; i++) {
            setIsTrait(i, false);
            useHeader(i,true);
        }
    }

    protected void setNumberOfTraits(int nTraits) {
        NumberOfTraits = nTraits;
        traitType = new String[NumberOfTraits];
        for (int i = 0; i < NumberOfTraits; i++) {
            traitType[i] = DATA;
        }
    }
}
