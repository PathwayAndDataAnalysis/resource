package org.panda.resource;

import org.biopax.paxtools.controller.Cloner;
import org.biopax.paxtools.controller.Completer;
import org.biopax.paxtools.controller.PathAccessor;
import org.biopax.paxtools.controller.SimpleEditorMap;
import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.BioPAXLevel;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.*;
import org.panda.utility.ContentSet;
import org.panda.utility.Kronometre;
import org.panda.utility.statistics.FDR;
import org.panda.utility.statistics.FishersExactTest;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class PCPathway extends FileServer
{
	private static PCPathway instance;
	private static final String FILE = "pcpathway.txt";

	private Map<String, Set<String>> gene2pathway;
	private Map<String, Set<String>> pathway2gene;
	private Map<String, String> pathway2name;
	private Map<String, String> pathway2resource;

	public static synchronized PCPathway get()
	{
		if (instance == null) instance = new PCPathway();
		return instance;
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{FILE};
	}

	@Override
	public boolean load() throws IOException
	{
		pathway2gene = new HashMap<>();
		gene2pathway = new HashMap<>();
		pathway2name = new HashMap<>();
		pathway2resource = new HashMap<>();

		Set<Set<String>> groups = new HashSet<>();

		Files.lines(Paths.get(locateInBase(getLocalFilenames()[0]))).map(line -> line.split("\t")).forEach(token ->
		{
			pathway2name.put(token[0], token[1]);
			pathway2resource.put(token[0], token[2]);

			if (token.length > 3)
			{
				Set<String> group = new ContentSet<>(Arrays.asList(token).subList(3, token.length));
				if (groups.contains(group)) return;
				groups.add(group);

				pathway2gene.put(token[0], group);

				for (int i = 2; i < token.length; i++)
				{
					if (!gene2pathway.containsKey(token[i])) gene2pathway.put(token[i], new HashSet<>());
					gene2pathway.get(token[i]).add(token[0]);
				}
			}
		});

		return true;
	}


	public Set<String> getPathways(String gene)
	{
		if (gene2pathway.containsKey(gene)) return gene2pathway.get(gene);
		else return Collections.emptySet();
	}

	public Set<String> getGenes(String pathwayID)
	{
		if (pathway2gene.containsKey(pathwayID)) return pathway2gene.get(pathwayID);
		else return Collections.emptySet();
	}

	public String getName(String id)
	{
		return pathway2name.get(id);
	}

	public String getResource(String id)
	{
		return pathway2resource.get(id);
	}

	public String getCoverageStr(String id, Set<String> genes)
	{
		if (!pathway2gene.containsKey(id)) return null;
		Set<String> set = new HashSet<String>(genes);
		set.retainAll(pathway2gene.get(id));
		return set.size() + "/" + pathway2gene.get(id).size();
	}

	/**
	 * Gets the enrichment pvals and pval limits.
	 * @param genes query
	 * @param background if there is any
	 * @return two maps, first is for pvals, second is for limits
	 */
	public Map<String, Double>[] getEnrichmentPvals(Collection<String> genes,
		Collection<String> background, int minMemberSize, int maxMemberSize)
	{
		if (background == null)
		{
			background = new HashSet<>(gene2pathway.keySet());
			if (!background.containsAll(genes))
			{
				Set<String> set = new HashSet<String>(genes);
				set.removeAll(background);

				genes = new HashSet<>(genes);
				genes.removeAll(set);
				System.out.println("Removed " + set.size() + " unknown genes: " + set);
				System.out.println("Using " + genes.size() + ": " + genes);
			}
		}

		if (!background.containsAll(genes)) throw new IllegalArgumentException(
			"Background genes have to contain all the selected genes.");

		Map<String, Integer> selectionCnt = count(genes);
		Map<String, Integer> backgroundCnt = count(background);

		Map<String, Double> mapP = new HashMap<>();
		Map<String, Double> mapL = new HashMap<>();

		for (String pathway : selectionCnt.keySet())
		{
			int size = pathway2gene.get(pathway).size();
			if (size < minMemberSize || size > maxMemberSize) continue;

			double pval = FishersExactTest.calcEnrichmentPval(background.size(),
				backgroundCnt.get(pathway), genes.size(), selectionCnt.get(pathway));

			double limit = FishersExactTest.calcEnrichmentPval(background.size(),
				backgroundCnt.get(pathway), genes.size(),
				Math.min(backgroundCnt.get(pathway), genes.size()));

			mapP.put(pathway, pval);
			mapL.put(pathway, limit);
		}

		return new Map[]{mapP, mapL};
	}

	private Map<String, Integer> count(Collection<String> genes)
	{
		Map<String, Integer> cnt = new HashMap<>();

		for (String pathway : pathway2gene.keySet())
		{
			Set<String> mems = new HashSet<>(pathway2gene.get(pathway));
			mems.retainAll(genes);
			if (!mems.isEmpty()) cnt.put(pathway, mems.size());
		}
		return cnt;
	}

	public List<String> getEnrichedPathways(Collection<String> genes,
		Collection<String> background, double fdrThr)
	{
		return getEnrichedPathways(genes, background, fdrThr, 3, 500);
	}

	public List<String> getEnrichedPathways(Collection<String> genes,
		Collection<String> background, double fdrThr, int minMemberSize, int maxMemberSize)
	{
		Map<String, Double>[] map = getEnrichmentPvals(genes, background, minMemberSize, maxMemberSize);
		if (fdrThr < 0)
		{
			fdrThr = FDR.decideBestFDR_BH(map[0], map[1]);
			System.out.println("fdrThr = " + fdrThr);
		}
		return FDR.select(map[0], map[1], fdrThr).stream().collect(Collectors.toList());
	}

	public Map<String, Double> getEnrichmentQvals(Collection<String> genes,
		Collection<String> background, int minMemberSize, int maxMemberSize)
	{
		Map<String, Double>[] map = getEnrichmentPvals(genes, background, minMemberSize, maxMemberSize);
		return FDR.getQVals(map[0], map[1]);
	}

	public void writeEnrichmentResults(Set<String> genes, int minMemberSize, int maxMemberSize, String filename)
		throws IOException
	{
		final Map<String, Double>[] pvals = getEnrichmentPvals(genes, null, minMemberSize, maxMemberSize);
		Map<String, Double> qvals = getEnrichmentQvals(genes, null, minMemberSize, maxMemberSize);

		List<String> ids = new ArrayList<>(qvals.keySet());
		Collections.sort(ids, (o1, o2) -> pvals[0].get(o1).compareTo(pvals[0].get(o2)));

		BufferedWriter writer = new BufferedWriter(new FileWriter(filename));

		writer.write("# Pathway name: Name of the pathway as given by the original pathway database.\n");
		writer.write("# P-value: Enrichment p-value of the pathway calculated by Fisher's exact test.\n");
		writer.write("# Q-value: Estimated FDR (false discovery rate) if this p-value is used as cutoff threshold.\n");
		writer.write("# Hit size: Number of query genes that overlaps with this pathway.\n");
		writer.write("# Pathway size: Number of genes in this pathway.\n");
		writer.write("# Genes contributed enrichment: Names of query genes that overlaps with this pathway.\n");
		writer.write("Pathway name\tResource\tPathway Commons ID\tP-value\tQ-value\tHit size\tPathway size\tGenes contributed enrichment");
		for (String id : ids)
		{
			if (pvals[0].get(id) > 0.05) break;

			Set<String> g = new HashSet<>(getGenes(id));
			int allSize = g.size();
			g.retainAll(genes);
			int hitSize = g.size();
			writer.write("\n" + getName(id) + "\t" + getResource(id) + "\t" + id + "\t" + pvals[0].get(id) + "\t" +
				qvals.get(id) + "\t" + hitSize + "\t" + allSize + "\t" + g);
		}
		writer.close();
	}

	/**
	 * Gets pathways sorted to their containment of the query genes. Does not control the density,
	 * but only the count, so it is likely that first pathways will be big ones.
	 */
	public TreeMap<String, Integer> getSortedPathways(Collection<String> genes)
	{
		final Map<String, Integer> map = new HashMap<String, Integer>();

		for (String gene : genes)
		{
			for (String pathway : getPathways(gene))
			{
				if (map.containsKey(pathway)) map.put(pathway, map.get(pathway) + 1);
				else map.put(pathway, 1);
			}
		}

		TreeMap<String, Integer> sorted = new TreeMap<String, Integer>((o1, o2) -> {
			return map.get(o2).compareTo(map.get(o1));
		});

		sorted.putAll(map);
		return sorted;
	}

	public static void main(String[] args) throws IOException
	{
//		String s = "PTEN, PIK3CA, ARID1A, PIK3R1, CTNNB1, TP53, KRAS, CTCF, FBXW7, LRP2, FGFR2, RYR1, TBL1XR1, MTOR, CACNA1A, PPP2R1A, PKN1, LYST, TRPM6, ERBB2, FN1, WDFY3, MYC, SPTB, DVL3, PRKCI, ECT2, ACTL6A, TBC1D31, IKBKB, PRKACA, DLG1, PTK2, THPO, DNM2, FOSL2, DSTYK, CCNE1, TNK2, EFNA1, PAK2, RASAL1, ARMC6, HGS, CDC37, TNFSF10, PPP1R1B, GRB2, PPP1CA";
//		List<String> select = getEnrichedPathways(Arrays.asList(s.split(", ")), null, 0.01);
//		for (String id : select)
//		{
//			System.out.println(id + "\t" + getName(id));
//		}

		prepareResource();
	}

	// Section: Preparing the pathways file from PC owl

	private static void prepareResource() throws IOException
	{
		Kronometre k = new Kronometre();
		Map<String, Set<String>> pathway2gene = new ConcurrentHashMap<>();
		Map<String, String> pathway2name = new ConcurrentHashMap<>();
		Map<String, String> pathway2resource = new ConcurrentHashMap<>();

		SimpleIOHandler h = new SimpleIOHandler();
		Model model = h.convertFromOWL(new FileInputStream("/home/babur/Documents/PC/PathwayCommons.8.Detailed.BIOPAX.owl"));

		System.out.println("Loaded BioPAX model");
		k.print();
		k.start();

		model.getObjects(Pathway.class).parallelStream().forEach(pathway ->
		{
			String id = pathway.getRDFId();
			String name = pathway.getDisplayName();

			if (name == null || name.isEmpty()) return;
			if (name.contains(" = ") || name.contains(" => ") || name.contains("[[") || name.contains("]]")) return;

			pathway2name.put(id, name);

			String resource = pathway.getDataSource().iterator().next().getName().iterator().next();
			pathway2resource.put(id, resource);

			Model m = excise(model, pathway);
			Set<String> syms = collectGeneSymbols(m);

//			Set<String> syms = collectSymbols(pathway, model);

			if (syms.size() < 3) return;

			if (!pathway2gene.containsKey(id)) pathway2gene.put(id, new HashSet<>());

			syms.forEach(pathway2gene.get(id)::add);
		});
		writeResources(pathway2gene, pathway2name, pathway2resource);
		k.print();
	}

	private static void writeResources(Map<String, Set<String>> pathway2gene, Map<String, String> pathway2name,
		Map<String, String> pathway2resource) throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter("../repo/resource-files/" + FILE));

		for (String id : pathway2gene.keySet())
		{
			if (pathway2gene.get(id).isEmpty()) continue;
			if (pathway2name.get(id).contains("$")) continue;

			writer.write(id + "\t" + pathway2name.get(id) + "\t" + pathway2resource.get(id));

			for (String sym : pathway2gene.get(id))
			{
				writer.write("\t" + sym);
			}
			writer.write("\n");
		}
		writer.close();
	}

	private static Model excise(Model model, Pathway pathway)
	{
		Completer c = new Completer(SimpleEditorMap.L3);

		Set<BioPAXElement> objects = c.complete(Collections.<BioPAXElement>singleton(pathway), model);

		Cloner cln = new Cloner(SimpleEditorMap.L3, BioPAXLevel.L3.getDefaultFactory());

		return cln.clone(model, objects);
	}

	private static Set<String> collectGeneSymbols(Model model)
	{
		Set<String> symbols = new HashSet<>();

		for (EntityReference er : model.getObjects(EntityReference.class))
		{
			for (Xref xref : er.getXref())
			{
				if (xref.getDb() == null) continue;
				if (xref.getDb().equalsIgnoreCase("HGNC SYMBOL"))
				{
					String s = HGNC.get().getSymbol(xref.getId());
					if (s != null) symbols.add(s);
				}
			}
		}
		return symbols;
	}

	private static Set<String> collectSymbols(Pathway pathway, Model model)
	{
		return new Completer(SimpleEditorMap.L3).complete(
			new HashSet<BioPAXElement>(new PathAccessor("Pathway/pathwayComponent*:Interaction").getValueFromBean(pathway)),
			model).stream()
				.filter(o -> o instanceof EntityReference).map(o -> (EntityReference) o)
				.map(XReferrable::getXref).flatMap(Collection::stream)
				.filter(xr -> xr.getDb() != null && xr.getDb().equalsIgnoreCase("hgnc symbol"))
				.map(xr -> HGNC.get().getSymbol(xr.getId())).filter(s -> s != null)
				.collect(Collectors.toSet());
	}
}
