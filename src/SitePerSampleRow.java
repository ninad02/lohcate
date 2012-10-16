//refName	coordinate	refBase	varBase	variantBase-N	variantBase-T	refEnv-N	refEnv-T	Q20_TotCov_N	Q20_TotCov_T	Q20_VarCov_N	Q20_VarCov_T	Q20_VariantRatio_N	Q20_VariantRatio_T	dbsnp	MutationType	Hugo_Symbol
//chr1	881627	G	A	A	A	AGGTCAGGGGTGT	AGGTCAGGGGTGT	75	72	75	71	1	0.986111	rs2272757,bycluster;byfrequency;by1000genomes,G|0.539;A|0.461,.	synonymous_SNV	NOC2L

public class SitePerSampleRow {

	/* We have a long member variable (0 - 63 bits) that stores several values encoded as bits:
	 * Bit  63: Kept as 0 so that it's not interpreted as a negative number
	 * Bits 58 - 62 ( 5 bits): chromosome number
	 * Bits 30 - 57 (28 bits): position on chromosome
	 * Bits 27 - 29 ( 3 bits): reference base
	 * Bits 24 - 26 ( 3 bits): variant base normal
	 * Bits 21 - 23 ( 3 bits): variant base tumor
	 * Bits 18 - 20 ( 3 bits): variant base (from new pipeline, sematics not clear)
	 * Bits 
	protected long mChrPos;
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
