package org.panda.resource;

import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * Provides genes in OncoKB.
 *
 * Update the second resource file by clicking the download link at https://www.oncokb.org/cancerGenes
 *
 * @author Ozgun Babur
 */
public class OncoKB extends FileServer
{
	private static OncoKB instance;

	private Map<String, Map<String, String[]>> variantData;
	private Set<String> cancerGenes;
	private Set<String> oncogenes;
	private Set<String> tumorSuppressors;

	public static synchronized OncoKB get()
	{
		if (instance == null) instance = new OncoKB();
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
		return new String[]{"oncoKB.txt", "OncoKB-cancerGeneList.tsv"};
	}

	@Override
	public String[] getDistantURLs()
	{
		return new String[]{"http://oncokb.org/api/v1/utils/allAnnotatedVariants.txt",
			GITHUB_REPO_BASE + "OncoKB-cancerGeneList.tsv"};
	}

	@Override
	public boolean load() throws IOException
	{
		variantData = new HashMap<>();

		getResourceAsStream(getLocalFilenames()[0], StandardCharsets.ISO_8859_1).skip(1).map(l -> l.split("\t"))
			.filter(t -> t.length > 1).forEach(t ->
		{
			if (!variantData.containsKey(t[0])) variantData.put(t[0], new HashMap<>());

			variantData.get(t[0]).put(t[1], ArrayUtil.getTail(t, 2));
		});

		String[] header = getResourceAsStream(getLocalFilenames()[1]).findFirst().get().split("\t");
		cancerGenes = new HashSet<>();
		oncogenes = new HashSet<>();
		tumorSuppressors = new HashSet<>();
		int oInd = ArrayUtil.indexOf(header, "Is Oncogene");
		int tInd = ArrayUtil.indexOf(header, "Is Tumor Suppressor Gene");

		getResourceAsStream(getLocalFilenames()[1]).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			cancerGenes.add(t[0]);
			if (t[oInd].equals("Yes")) oncogenes.add(t[0]);
			if (t[tInd].equals("Yes")) tumorSuppressors.add(t[0]);
		});

		return true;
	}


	public static void main(String[] args) throws IOException
	{
		CollectionUtil.printVennCounts(get().oncogenes, get().tumorSuppressors);
//		writeHighlightFile("/Users/ozgun/Documents/Analyses/CPTAC-LSCC-3.2/oncokb.highlight");
	}

	//----SMMART related-----
	public static void annotateSMMARTMutations(String[] args) throws IOException
	{
		Map<String, Set<String>> map = new HashMap<>();

		Files.lines(Paths.get("/home/babur/Documents/Analyses/SMMART/Patient1-revisit/variantData/met1-muts.maf"))
			.map(l -> l.split("\t")).filter(t -> t.length > 41).filter(t -> !t[41].isEmpty()).forEach(t ->
		{
			if (!map.containsKey(t[0])) map.put(t[0], new HashSet<>());
			map.get(t[0]).add(t[41].substring(2));
		});

		get().printMatchReport(map);
	}

	public void printMatchReport(Map<String, Set<String>> muts)
	{
		muts.keySet().stream().sorted().forEach(gene ->
		{
			if (variantData.containsKey(gene))
			{
				for (String mut : muts.get(gene))
				{
					if (variantData.get(gene).containsKey(mut))
					{
						System.out.println(gene + "\t" + mut + "\t" +  Arrays.toString(variantData.get(gene).get(mut)));
					}
					else if ((mut.endsWith("*") || mut.contains("fs")) &&
						variantData.get(gene).containsKey("Truncating Mutations"))
					{
						System.out.println(gene + "\t" + mut + " (trunc)"  + "\t" +  Arrays.toString(variantData.get(gene).get("Truncating Mutations")));
					}
					else
					{
						System.out.println(gene + "\t" + mut);
					}
				}
			}
		});
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
