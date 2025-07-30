package org.panda.resource;

import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashSet;
import java.util.Set;

/**
 * To parse COSMIC hallmarks file
 * @author Ozgun Babur
 */
public class COSMICLungHallmarkCancerGenes extends FileServer
{
	private static COSMICLungHallmarkCancerGenes instance;

	private Set<String> cancerGenes;
	private Set<String> oncogenes;
	private Set<String> tumorSuppressors;

	public static synchronized COSMICLungHallmarkCancerGenes get()
	{
		if (instance == null) instance = new COSMICLungHallmarkCancerGenes();
		return instance;
	}

	public Set<String> getAllSymbols()
	{
		return cancerGenes;
	}

	public Set<String> getCancerGenes()
	{
		return cancerGenes;
	}

	public boolean isCancerGene(String sym)
	{
		return cancerGenes.contains(sym);
	}

	public boolean isOncogene(String sym)
	{
		return oncogenes.contains(sym);
	}

	public boolean isOncogeneOnly(String sym)
	{
		return isOncogene(sym) && !isTumorSuppressor(sym);
	}

	public boolean isTumorSuppressor(String sym)
	{
		return tumorSuppressors.contains(sym);
	}

	public boolean isTumorSuppressorOnly(String sym)
	{
		return isTumorSuppressor(sym) && !isOncogene(sym);
	}

	public Set<String> getOncogenes()
	{
		return this.oncogenes;
	}

	public Set<String> getTumorSuppressors()
	{
		return this.tumorSuppressors;
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"Cosmic_CancerGeneCensusHallmarksOfCancer_v102_GRCh38.tsv"};
	}

	@Override
	public boolean load() throws IOException
	{
		cancerGenes = new HashSet<>();
		oncogenes = new HashSet<>();
		tumorSuppressors = new HashSet<>();

		String[] header = getResourceAsStream(getLocalFilenames()[0]).findFirst().get().split("\t");
		int symInd = ArrayUtil.indexOf(header, "GENE_SYMBOL");
		int cellTypeInd = ArrayUtil.indexOf(header, "CELL_TYPE");
		int hallmarkInd = ArrayUtil.indexOf(header, "HALLMARK");
		int impactInd = ArrayUtil.indexOf(header, "IMPACT");

		getResourceAsStream(getLocalFilenames()[0]).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String symbol = t[symInd];
			String cellType = t[cellTypeInd];

			if (cellType.contains("lung") || cellType.contains("NSCLC"))
			{
				String hallmark = t[hallmarkInd];
				if (hallmark.equals("role in cancer"))
				{
					String impact = t[impactInd];
					if (impact.contains("oncogene") || impact.contains("TSG"))
					{
						cancerGenes.add(symbol);

						if (impact.contains("oncogene")) oncogenes.add(symbol);
						if (impact.contains("TSG")) tumorSuppressors.add(symbol);
					}
				}
			}
		});

		return true;
	}


	public static void main(String[] args) throws IOException
	{
		CollectionUtil.printVennSets(get().oncogenes, get().tumorSuppressors, GeminiLungCancerGenes.get().getOncogenes(), GeminiLungCancerGenes.get().getTumorSuppressors());
//		writeHighlightFile("/Users/ozgun/Documents/Analyses/CPTAC-LSCC-3.2/oncokb.highlight");
	}

	private static void writeHighlightFile(String file) throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(file));

		for (String sym : get().getCancerGenes())
		{
			writer.write("node\t" + sym + "\n");
		}

		writer.close();
	}
}
