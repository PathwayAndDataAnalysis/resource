package org.panda.resource;

import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.SIFFileUtil;
import org.panda.utility.TermCounter;
import org.panda.utility.statistics.FishersExactTest;
import org.panda.utility.statistics.GeneSetEnrichment;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Provides GO ontology graph.
 *
 * @author Ozgun Babur
 */
public class GO extends FileServer
{
	private Map<String, String> idToName;
	private Map<String, Set<String>> geneToGO;
	private Map<String, Set<String>> goToGene;
	private Map<String, Set<String>> isAMap;

	private static GO instance;

	public String getNameOfTerm(String goID)
	{
		return idToName.get(goID);
	}

	public Set<String> getGOIDs(String gene)
	{
		if (geneToGO.containsKey(gene)) return geneToGO.get(gene);
		return Collections.emptySet();
	}

	public Set<String> getCommonGOIDs(String... gene)
	{
		if (gene == null) return Collections.emptySet();
		if (gene.length == 1) return getGOIDs(gene[0]);

		Set<String> ids = getGOIDs(gene[0]);
		for (int i = 1; i < gene.length; i++)
		{
			ids.retainAll(getGOIDs(gene[i]));
		}
		return ids;
	}

	public Set<String> getAllGenes()
	{
		return geneToGO.keySet();
	}

	public Set<String> getGenes(Set<String> terms)
	{
		Set<String> genes = new HashSet<>();
		for (String term : terms)
		{
			genes.addAll(getGenes(term));
		}
		return genes;
	}

	public Set<String> getGenes(String term)
	{
		if (goToGene.containsKey(term)) return goToGene.get(term);
		return Collections.emptySet();
	}

	public Set<String> getParentTerms(String term)
	{
		if (isAMap.containsKey(term)) return isAMap.get(term);
		return Collections.emptySet();
	}

	public void printAssociatedTerms(Set<String> genes)
	{
		printAssociatedTerms(genes, Collections.emptySet(), Collections.emptySet());
	}

	public void printTermsOfKeyword(String keyword)
	{
		Map<String, String> map = getTermsContaining(keyword);
		map.forEach((id,  name) -> System.out.println(id + "\t" + name));
	}

	public void printAssociatedTerms(Set<String> genes, Set<String> ignoreTerms, Set<String> ignoreGenesOfTerms)
	{
		Set<String> termSet = genes.stream().map(this::getGOIDs).flatMap(Collection::stream).collect(Collectors.toSet());
		termSet.removeAll(ignoreTerms);

		genes = new HashSet<>(genes);
		genes.removeAll(getGenes(ignoreGenesOfTerms));

		Map<String, Set<String>> map = new HashMap<>();
		for (String term : termSet)
		{
			Set<String> set = new HashSet<>(getGenes(term));
			set.retainAll(genes);
			assert !set.isEmpty();
			map.put(term, set);
		}
		termSet.stream().sorted((t1, t2) -> new Integer(map.get(t2).size()).compareTo(map.get(t1).size()))
			.forEach(term -> System.out.println(term + "\t" + getNameOfTerm(term) + "\t" + map.get(term)));
	}

	public void printAssociatedTerms(String gene, Set<String> ignore)
	{
		System.out.println("\ngene = " + gene);

		if (geneToGO.containsKey(gene))
		{
			for (String term : geneToGO.get(gene))
			{
				if (ignore == null || !ignore.contains(term))
				{
					String name = idToName.get(term);

					if (!name.startsWith("positive regulation ") && !name.startsWith("negative regulation "))
					{
						System.out.println(term + "\t" + name);
					}
				}
			}
		}
	}

	public void printAssociatedCommonTerms(String... gene)
	{
		Set<String> ids = getCommonGOIDs(gene);

		for (String id : ids)
		{
			String name = idToName.get(id);

			if (!name.startsWith("positive regulation ") && !name.startsWith("negative regulation "))
			{
				System.out.println(id + "\t" + name);
			}
		}
	}

	public Set<String> getGenesOfTerm(String goID)
	{
		return goToGene.get(goID);
	}

	public Map<String, String> getTermsContaining(String query)
	{
		return idToName.keySet().stream().filter(id -> idToName.get(id).toLowerCase().contains(query.toLowerCase()))
			.collect(Collectors.toMap(id -> id, id -> idToName.get(id)));
	}

	public Set<String> getGenesContainingKeywordInTermNames(String keyword)
	{
		return geneToGO.keySet().stream().filter(g -> geneToGO.get(g).stream().filter(idToName::containsKey).map(idToName::get)
			.anyMatch(name -> name.contains(keyword))).collect(Collectors.toSet());
	}

	public Set<String> getGenesContainingKeywordInTermNames(Set<String> keywords)
	{
		return getGenesContainingKeywordInTermNames(keywords, Collections.emptySet());
	}

	/**
	 * Finds genes that has at least one GO term that contains one of the keywordsToKeep, but that same GO term does not
	 * contain any keywordsToSkip.
	 */
	public Set<String> getGenesContainingKeywordInTermNames(Set<String> keywordsToKeep, Set<String> keywordsToSkip)
	{
		return geneToGO.keySet().stream().filter(g -> geneToGO.get(g).stream().map(idToName::get)
			.anyMatch(name -> keywordsToKeep.stream().anyMatch(name::contains) &&
				keywordsToSkip.stream().noneMatch(name::contains)))
			.collect(Collectors.toSet());
	}

	public void writeTermsContainingKeywordForGivenGenes(Set<String> genes, Set<String> keywordsToKeep, Set<String> keywordsToSkip, String outFile)
	{
		Map<String, Set<String>> goToSubset = new HashMap<>();

		idToName.keySet().stream().filter(id -> goToGene.containsKey(id) && !goToGene.get(id).isEmpty())
			.filter(id -> contains(idToName.get(id), keywordsToKeep) && !contains(idToName.get(id), keywordsToSkip))
			.forEach(id ->
		{
			Set<String> hit = CollectionUtil.getIntersection(genes, goToGene.get(id));
			if (!hit.isEmpty())
			{
				goToSubset.put(id, hit);
			}
		});

		List<String> ids = new ArrayList<>(goToSubset.keySet());
		ids.sort((o1, o2) -> Integer.compare(goToSubset.get(o2).size(), goToSubset.get(o1).size()));

		BufferedWriter writer = FileUtil.newBufferedWriter(outFile);
		FileUtil.write("ID\tTerm\tGenes", writer);

		for (String id : ids)
		{
			FileUtil.lnwrite(id + "\t" + idToName.get(id) + "\t", writer);
			goToSubset.get(id).stream().sorted().forEach(gene -> FileUtil.write(gene + " ", writer));
		}

		FileUtil.closeWriter(writer);
	}

	public boolean contains(String s, Set<String> q)
	{
		for (String query : q)
		{
			if (s.contains(query)) return true;
		}
		return false;
	}

	/**
	 * Returns Fisher's exact test p-values for enrichment of each GO term. Does not correct for multiple hypothesis
	 * testing.
	 *
	 * @param background pass null if the background is just all possible genes
	 */
	public Map<String, Double> calculateEnrichment(Collection<String> queryGenes, Collection<String> background,
		int minNumOfMembersInAGroup, int maxNumOfMembersInAGroup)
	{
		return GeneSetEnrichment.calculateEnrichment(queryGenes, background, minNumOfMembersInAGroup,
			maxNumOfMembersInAGroup, goToGene);
	}


	public static GO get()
	{
		if (instance == null) instance = new GO();
		return instance;
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"GO-ontology.txt", "GO-gene-association.txt"};
	}

	@Override
	public String[] getDistantURLs()
	{
		return new String[]{"http://purl.obolibrary.org/obo/go/go-basic.obo",
			"http://geneontology.org/gene-associations/goa_human.gaf.gz"};
	}

	@Override
	public boolean load() throws IOException
	{
		idToName = new HashMap<>();
		isAMap = new HashMap<>();
		Scanner sc = new Scanner(new File(ResourceDirectory.get() + File.separator + getLocalFilenames()[0]));
		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			if (line.startsWith("id: "))
			{
				String term = line.substring(4);
				idToName.put(term, sc.nextLine().substring(6));

				while (sc.hasNextLine() && !line.startsWith("["))
				{
					line = sc.nextLine();
					if (line.startsWith("is_a: "))
					{
						String parent = line.substring(6);
						parent = parent.substring(0, parent.indexOf(" "));
						if (!isAMap.containsKey(term)) isAMap.put(term, new HashSet<>());
						isAMap.get(term).add(parent);
					}
				}
			}
		}

		geneToGO = new HashMap<>();
		goToGene = new HashMap<>();

		getResourceAsStream(getLocalFilenames()[1]).filter(l -> !l.startsWith("!")).map(l -> l.split("\t")).forEach(t -> {
			String gene = t[2];
			String rel = t[3];
			String go = t[4];

			if (rel.contains("NOT")) return;

			if (!geneToGO.containsKey(gene)) geneToGO.put(gene, new HashSet<>());
			geneToGO.get(gene).add(go);
			if (!goToGene.containsKey(go)) goToGene.put(go, new HashSet<>());
			goToGene.get(go).add(gene);
		});

		for (String term : new HashSet<>(goToGene.keySet()))
		{
			Set<String> genes = goToGene.get(term);

			for (String parent : getParentTerms(term))
			{
				addParentsRecursive(parent, genes);
			}
		}

		return true;
	}

	private void addParentsRecursive(String parent, Set<String> genes)
	{
		if (!goToGene.containsKey(parent)) goToGene.put(parent, new HashSet<>());
		goToGene.get(parent).addAll(genes);

		if (genes.contains("SCD") && parent.equals("GO:0046390"))
		{
			System.out.println();
		}



		for (String gene : genes)
		{
			geneToGO.get(gene).add(parent);
		}

		for (String grand : getParentTerms(parent))
		{
			addParentsRecursive(grand, genes);
		}
	}

	// Section: Temporary methods

	public static void main(String[] args) throws IOException
	{
//		List<String> genes = readPanCanGenes();
//		Set<String> ignore = readIgnoreList();

//		get().printAssociatedTerms(new HashSet<>(genes));

//		genes.forEach(g -> get().printAssociatedTerms(g, ignore));
//		get().printAssociatedTerms("KDM5A", Collections.emptySet());

//		checkInterestTermAssociationRate();

//		String go = "GO:0046390";
//		System.out.println(go + " = " + get().getNameOfTerm(go));

//		get().printAssociatedCommonTerms("SCD");

//		Set<String> genes1 = get().getGenesOfTerm("GO:0006629");
//		Set<String> genes2 = get().getGenesOfTerm("GO:0008203");
//		CollectionUtil.printVennCounts(genes1, genes2);
//		genes2.stream().sorted().forEach(System.out::println);

		get().printTermsOfKeyword("calcium");
//		get().printTermsOfKeyword("histone acetyltransferase activity");
		get().getGenesContainingKeywordInTermNames("calcium").stream().sorted().forEach(System.out::println);
//		get().getGenesContainingKeywordInTermNames("regulation of histone deacetylase activity").stream().sorted().forEach(System.out::println);

//		get().getGenesContainingKeywordInTermNames(new HashSet<>(Arrays.asList("histone acetyltransferase activity", "histone deacetylase activity")), new HashSet<>(Arrays.asList("regulation of"))).forEach(System.out::println);

//		Set<String> genes = SIFFileUtil.getGenesInSIFFile("/home/ozgunbabur/Downloads/temp/causative-network-only.sif");
//		get().writeTermsContainingKeywordForGivenGenes(genes, new HashSet<>(Arrays.asList("calcium", "viral", "virus")), Collections.emptySet(), "/home/ozgunbabur/Downloads/temp/go.txt");

//		get().printAverageSetSize();
	}

	private static List<String> readPanCanGenes() throws IOException
	{
		return readFirstColumn("/home/babur/Documents/PanCan/pancan.txt", 20);
	}
	private static Set<String> readIgnoreList() throws IOException
	{
		return new HashSet<>(readFirstColumn("/home/babur/Documents/PanCan/GO-term-ignore-list.txt", Integer.MAX_VALUE));
	}

	private static List<String> readFirstColumn(String file, int lineLimit) throws IOException
	{
		return Files.lines(Paths.get(file)).skip(1).limit(lineLimit)
			.map(l -> l.split("\t")[0]).collect(Collectors.toList());
	}

	private static void checkInterestTermAssociationRate() throws IOException
	{
		Set<String> interest = Files.lines(Paths.get("/home/babur/Documents/PanCan/GO-term-interest.txt"))
			.filter(l -> !l.startsWith("#")).collect(Collectors.toSet());

		System.out.println("Total genes = " + get().geneToGO.keySet().size());

		for (String gene : get().geneToGO.keySet())
		{
			for (String id : get().geneToGO.get(gene))
			{
				if (get().idToName.get(id) == null) System.err.println("Hell!");
			}
		}

		Set<String> match = new HashSet<>();
		for (String s : interest)
		{
			Set<String> genes = get().getGenesContainingKeywordInTermNames(s);
			int cnt = genes.size();
			System.out.println(cnt + "\t" + s);
			match.addAll(genes);
		}
		System.out.println();
		System.out.println("match.size() = " + match.size());
	}

	public void printAverageSetSize()
	{
		int[] cc = new int[]{0, 0};
		TermCounter tc = new TermCounter();
		goToGene.forEach((s, strings) ->
		{
			if (!strings.isEmpty())
			{
				cc[0]++;
				cc[1] += strings.size();

				for (String gene : strings)
				{
					tc.addTerm(gene);
				}
			}
		});

		System.out.println(cc[1] / (double) cc[0]);
		System.out.println("cc = " + Arrays.toString(cc));
		System.out.println("tc.getMeanCountPerTerm() = " + tc.getMeanCountPerTerm());
	}
}
