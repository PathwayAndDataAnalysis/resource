package org.panda.resource;

import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Stream;

/**
 * This class provides a mapping between HGNC IDs and approved gene symbols.
 *
 * @author Ozgun Babur
 */
public class MGIVertebrateHomology extends FileServer
{
	private static MGIVertebrateHomology instance;

	private Map<Organism, Map<String, Set<String>>> other2human;
	private Map<Organism, Map<String, Set<String>>> human2other;

	public static synchronized MGIVertebrateHomology get()
	{
		if (instance == null) instance = new MGIVertebrateHomology();
		return instance;
	}

	public Set<String> getCorrespondingHumanSymbols(String sym, Organism org)
	{
		return other2human.get(org).getOrDefault(sym, Collections.emptySet());
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"mgi-vertebrate-homology.txt"};
	}

	@Override
	public String[] getDistantURLs()
	{
		return new String[]{"http://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt"};
	}

	@Override
	public boolean load() throws IOException
	{
		human2other = new HashMap<>();
		other2human = new HashMap<>();

		Scanner scanner = new Scanner(new File(locateInBase(getLocalFilenames()[0])));
		int group = 0;

		scanner.nextLine();

		Map<Organism, Set<String>> symbolMap = new HashMap<>();

		while (scanner.hasNextLine())
		{
			String[] t = scanner.nextLine().split("\t");

			int g = Integer.valueOf(t[0]);
			if (g != group)
			{
				if (symbolMap.size() > 1 && symbolMap.containsKey(Organism.HUMAN))
				{
					Set<String> hSet = symbolMap.get(Organism.HUMAN);

					for (Organism org : symbolMap.keySet())
					{
						if (org != Organism.HUMAN)
						{
							Set<String> oSet = symbolMap.get(org);
							if (!human2other.containsKey(org)) human2other.put(org, new HashMap<>());
							if (!other2human.containsKey(org)) other2human.put(org, new HashMap<>());

							for (String h : hSet)
							{
								for (String m : oSet)
								{
									if (!human2other.get(org).containsKey(h))
										human2other.get(org).put(h, new HashSet<>());
									if (!other2human.get(org).containsKey(m))
										other2human.get(org).put(m, new HashSet<>());
									human2other.get(org).get(h).add(m);
									other2human.get(org).get(m).add(h);
								}
							}
						}
					}
				}

				symbolMap.clear();
				group = g;
			}

			if (t[1].equals("human"))
			{
				if (!symbolMap.containsKey(Organism.HUMAN)) symbolMap.put(Organism.HUMAN, new HashSet<>());
				symbolMap.get(Organism.HUMAN).add(t[3]);
			}
			else
			{
				Organism org = Organism.get(t[1]);
				assert org != null;

				if (!symbolMap.containsKey(org)) symbolMap.put(org, new HashSet<>());
				symbolMap.get(org).add(t[3]);
			}
		}

		// no need to handle the last row since we know it has no human corresponding gene
		return true;
	}

	public enum Organism
	{
		HUMAN("human"),
		MOUSE("mouse, laboratory"),
		RAT("rat"),
		CHIMPANZEE("chimpanzee"),
		MACAQUE("macaque, rhesus"),
		DOG("dog, domestic"),
		CATTLE("cattle"),
		FROG("frog, western clawed"),
		ZEBRAFISH("zebrafish"),
		CHICKEN("chicken");

		String name;

		Organism(String name)
		{
			this.name = name;
		}

		static Organism get(String name)
		{
			for (Organism org : values())
			{
				if (org.name.equals(name)) return org;
			}
			return null;
		}
	}


	public static void main(String[] args)
	{
		Set<String> terms = FileUtil.getTermsInTabDelimitedColumn("/home/ozgunbabur/Data/Linden/LiverRNA/RNA seq/DESeq2/5_stz_veh_hfd_sed-stz_sglt2i_hfd_sed_31 removed/DESeq2_output/DESeq2_stats_allGenes_CookTested.txt", 0, 1);
		Map<String, String> map1 = new HashMap<>();
		terms.stream().map(s -> s.replaceAll("\"", "")).forEach(s ->
		{
			Set<String> set = get().getCorrespondingHumanSymbols(s, Organism.RAT);
			if (set.size() == 1) map1.put(s, set.iterator().next());
			else if (set.contains(s.toUpperCase())) map1.put(s, s.toUpperCase());
		});
		System.out.println("map1 = " + map1.size());
		Map<String, Integer> cMap = new HashMap<>();
		map1.forEach((s1, s2) -> cMap.put(s2, cMap.getOrDefault(s2, 0) + 1));
		long count = map1.keySet().stream().filter(s -> cMap.get(map1.get(s)) == 1).count();
		System.out.println("count = " + count);
	}
}
