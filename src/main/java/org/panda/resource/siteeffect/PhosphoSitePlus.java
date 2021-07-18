package org.panda.resource.siteeffect;

import org.panda.resource.HGNC;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Stream;

/**
 * Provides phospho-site effects from PhosphoSitePlus.
 *
 * @author Ozgun Babur
 */
public class PhosphoSitePlus extends SiteEffectServer
{
	private static PhosphoSitePlus instance;

	private Map<Feature, Map<String, Map<String, String>>> actualMap;

	public static synchronized PhosphoSitePlus get()
	{
		if (instance == null) instance = new PhosphoSitePlus();
		return instance;
	}

	public static synchronized void initSingletonWith(Stream<String> resourceStream)
	{
		instance = new PhosphoSitePlus()
		{
			@Override
			public synchronized boolean init() throws IOException
			{
				return true;
			}
		};

		instance.loadFrom(resourceStream);
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"psp-regulatory-sites", "manually-curated-sites.txt"};
	}

	protected void printSites(String gene)
	{
		System.out.println("Gene: " + gene);
		for (Feature mod : Feature.values())
		{
			System.out.println("mod = " + mod);
			if (typeMap.get(mod).containsKey(gene))
			{
				for (String site : sortSites(typeMap.get(mod).get(gene).keySet()))
				{
					Integer sign = typeMap.get(mod).get(gene).get(site);
					System.out.print("\tsite: " + site + "\t" + (sign == 1 ? "activating" : sign == -1 ? "inhibiting" : "complex"));
					System.out.println("\t(" + actualMap.get(mod).get(gene).get(site) + ")");
				}
			}
			else
			{
				System.out.println("Not found.");
			}
		}
	}

	@Override
	public boolean load() throws IOException
	{
		typeMap = new HashMap<>();
		actualMap = new HashMap<>();

		Files.lines(Paths.get(locateInBase(getLocalFilenames()[0]))/*, Charset.forName("windows-31j")*/).skip(4)
			.map(line -> line.split("\t"))
			.filter(token -> token.length >= 13 && token[6].equals("human") &&
				HGNC.get().getSymbol(token[0]) != null)
			.forEach(token ->
			{
				String gene = HGNC.get().getSymbol(token[0]);
				String site = token[7];

				String modStr = site.substring(site.lastIndexOf("-") + 1);
				site = site.substring(0, site.lastIndexOf("-"));

				Feature mod = null;

				switch (modStr)
				{
					case "p":
						mod = Feature.PHOSPHORYLATION;
						break;
					case "ac":
						mod = Feature.ACETYLATION;
						break;
					case "me":
					case "m1":
					case "m2":
						mod = Feature.METHYLATION;
						break;
					case "ub":
						mod = Feature.UBIQUITINATION;
						break;
				}
				if (mod == null) return;

				if (!typeMap.containsKey(mod))
				{
					typeMap.put(mod, new HashMap<>());
					actualMap.put(mod, new HashMap<>());
				}

				if (!typeMap.get(mod).containsKey(gene)) typeMap.get(mod).put(gene, new HashMap<>());
				if (!actualMap.get(mod).containsKey(gene)) actualMap.get(mod).put(gene, new HashMap<>());


				if (typeMap.get(mod).get(gene).containsKey(site)) return;

				String function = token[11];
				actualMap.get(mod).get(gene).put(site, function);

				boolean actWord = false;
				boolean inhWord = false;

				if ((function.contains("induced") &&
					!function.contains("receptor desensitization, induced")))
				{
					actWord = true;
				}
				if ((function.contains("inhibited") ||
					function.contains("receptor desensitization, induced")))
				{
					inhWord = true;
				}

				if (actWord == inhWord)
				{
					if (function.contains("stabilization"))
					{
						actWord = true;
					}
					if (function.contains("degradation"))
					{
						inhWord = true;
					}
				}

				if (actWord == inhWord)
				{
					typeMap.get(mod).get(gene).put(site, 0);
				} else
				{
					typeMap.get(mod).get(gene).put(site, actWord ? 1 : -1);
				}
		});

		Files.lines(Paths.get(locateInBase(getLocalFilenames()[1]))).filter(l -> !l.startsWith("#"))
			.map(line -> line.split("\\s+")).filter(token -> token.length > 2).forEach(token ->
		{
			String gene = token[0];

			if (!typeMap.get(Feature.PHOSPHORYLATION).containsKey(gene)) typeMap.get(Feature.PHOSPHORYLATION).put(gene, new HashMap<>());
			if (!actualMap.get(Feature.PHOSPHORYLATION).containsKey(gene)) actualMap.get(Feature.PHOSPHORYLATION).put(gene, new HashMap<>());

			String site = token[1];
			int sign = Integer.parseInt(token[2]);

			typeMap.get(Feature.PHOSPHORYLATION).get(gene).put(site, sign);
			actualMap.get(Feature.PHOSPHORYLATION).get(gene).put(site, "manual curation");
		});

		return true;
	}

	public void loadFrom(Stream<String> resourceStream)
	{
		typeMap = new HashMap<>();

		resourceStream.map(l -> l.split("\t")).forEach(t ->
		{
			Feature mod = Feature.valueOf(t[2]);
			if (typeMap.containsKey(mod)) typeMap.put(mod, new HashMap<>());

			if (!typeMap.get(mod).containsKey(t[0])) typeMap.get(mod).put(t[0], new HashMap<>());
			typeMap.get(mod).get(t[0]).put(t[1], Integer.valueOf(t[3]));
		});
	}

	public static void main(String[] args)
	{
		PhosphoSitePlus psp = new PhosphoSitePlus();
//		List<String> list = getGenesWithMostSites();
//		for (int i = 0; i < 10; i++)
//		{
//			printSites(list.get(i));
//		}
		psp.printSites("FOXM1");
//		printUniqueAA();

//		List<Integer> dists = new ArrayList<>();
//		for (String gene : typeMap.keySet())
//		{
//			Map<String, Integer> sites = typeMap.get(gene);
//			int min = Integer.MAX_VALUE;
//
//			for (String s1 : sites.keySet())
//			{
//				for (String s2 : sites.keySet())
//				{
//					if (sites.get(s1) * sites.get(s2) == -1)
//					{
//						int dif = Math.abs(Integer.parseInt(s1.substring(1)) - Integer.parseInt(s2.substring(1)));
//						if (dif < min) min = dif;
//					}
//				}
//			}
//
//			if (min < Integer.MAX_VALUE) dists.add(min);
//			if (min < 10)
//			{
//				System.out.println("\n" + gene + "\t" + min);
//				printSites(gene);
//			}
//		}
//
//		Histogram h = new Histogram(10);
//		h.setBorderAtZero(true);
//		for (Integer dist : dists)
//		{
//			h.count(dist);
//		}
//		h.print();
	}

}
