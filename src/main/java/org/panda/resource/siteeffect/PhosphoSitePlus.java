package org.panda.resource.siteeffect;

import org.panda.resource.HGNC;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * Provides phospho-site effects from PhosphoSitePlus.
 *
 * @author Ozgun Babur
 */
public class PhosphoSitePlus extends SiteEffectServer
{
	private static PhosphoSitePlus instance;

	Map<String, Map<String, String>> actualMap;

	public static synchronized PhosphoSitePlus get()
	{
		if (instance == null) instance = new PhosphoSitePlus();
		return instance;
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"Regulatory_sites", "manually-curated-sites.txt"};
	}

	protected void printSites(String gene)
	{
		System.out.println("Gene: " + gene);
		if (typeMap.containsKey(gene))
		{
			for (String site : sortSites(typeMap.get(gene).keySet()))
			{
				Integer sign = typeMap.get(gene).get(site);
				System.out.print("\tsite: " + site + "\t" + (sign == 1 ? "activating" : sign == -1 ? "inhibiting" : "complex"));
				System.out.println("\t(" + actualMap.get(gene).get(site) + ")");
			}
		}
		else
		{
			System.out.println("Not found.");
		}
	}

	@Override
	public boolean load() throws IOException
	{
		typeMap = new HashMap<>();
		actualMap = new HashMap<>();

		Files.lines(Paths.get(locateInBase(getLocalFilenames()[1]))).filter(l -> !l.startsWith("#"))
			.map(line -> line.split("\\s+")).filter(token -> token.length > 2).forEach(token ->
		{
			String gene = token[0];

			if (!typeMap.containsKey(gene)) typeMap.put(gene, new HashMap<>());
			if (!actualMap.containsKey(gene)) actualMap.put(gene, new HashMap<>());

			String site = token[1];
			int sign = Integer.parseInt(token[2]);

			typeMap.get(gene).put(site, sign);
			actualMap.get(gene).put(site, "manual curation");
		});

		Files.lines(Paths.get(locateInBase(getLocalFilenames()[0]))/*, Charset.forName("windows-31j")*/).skip(4)
			.map(line -> line.split("\t"))
			.filter(token -> token.length >= 13 && token[6].equals("human") &&
				token[8].equals("PHOSPHORYLATION") && HGNC.get().getSymbol(token[4]) != null)
			.forEach(token -> {
				String gene = HGNC.get().getSymbol(token[4]);
				if (!typeMap.containsKey(gene)) typeMap.put(gene, new HashMap<>());
				if (!actualMap.containsKey(gene)) actualMap.put(gene, new HashMap<>());

				String site = token[7];

				if (typeMap.get(gene).containsKey(site)) return;

				actualMap.get(gene).put(site, token[12]);

				boolean actWord = false;
				boolean inhWord = false;

				if ((token[12].contains("induced") &&
					!token[12].contains("receptor desensitization, induced")))
				{
					actWord = true;
				}
				if ((token[12].contains("inhibited") ||
					token[12].contains("receptor desensitization, induced")))
				{
					inhWord = true;
				}

				if (actWord == inhWord)
				{
					if (token[12].contains("stabilization"))
					{
						actWord = true;
					}
					if (token[12].contains("degradation"))
					{
						inhWord = true;
					}
				}

				if (actWord == inhWord)
				{
					typeMap.get(gene).put(site, 0);
				} else
				{
					typeMap.get(gene).put(site, actWord ? 1 : -1);
				}
		});

		return true;
	}

	public static void main(String[] args)
	{
		PhosphoSitePlus psp = new PhosphoSitePlus();
//		List<String> list = getGenesWithMostSites();
//		for (int i = 0; i < 10; i++)
//		{
//			printSites(list.get(i));
//		}
		psp.printSites("FGR");
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
