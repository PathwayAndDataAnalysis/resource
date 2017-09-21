package org.panda.resource.siteeffect;

import org.panda.resource.FileServer;
import org.panda.resource.tcga.ProteomicsFileRow;
import org.panda.utility.TermCounter;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * Provides phospho-site effects from PhosphoSitePlus.
 *
 * @author Ozgun Babur
 */
public abstract class SiteEffectServer extends FileServer
{
	Map<String, Map<String, Integer>> typeMap;

	public Integer getEffect(String gene, String site)
	{
		if (typeMap.containsKey(gene))
		{
			return typeMap.get(gene).get(site);
		}
		return null;
	}

	public Integer getClosestEffect(String gene, String site, int distanceThreshold)
	{
		if (typeMap.containsKey(gene))
		{
			int s0 = Integer.parseInt(site.substring(1));

			Integer effect = null;
			int closestDist = Integer.MAX_VALUE;

			for (String ss : typeMap.get(gene).keySet())
			{
				Integer eff = typeMap.get(gene).get(ss);

				int s1 = Integer.parseInt(ss.substring(1));

				int dist = Math.abs(s0 - s1);
				if (dist == closestDist && !eff.equals(effect))
				{
					effect = 0;
				}
				else if (dist < closestDist)
				{
					effect = eff;
					closestDist = dist;
				}
			}
			if (closestDist <= distanceThreshold) return effect;
		}
		return null;
	}

	protected List<String> sortSites(Set<String> sites)
	{
		List<String> list = new ArrayList<>(sites);
		Collections.sort(list, (o1, o2) -> {
			try
			{
				return new Integer(o1.substring(1)).compareTo(new Integer(o2.substring(1)));
			}
			catch (NumberFormatException e)
			{
				return 0;
			}
		});
		return list;
	}

	protected List<String> getGenesWithMostSites()
	{
		List<String> genes = new ArrayList<>(typeMap.keySet());
		Collections.sort(genes, (o1, o2) -> new Integer(typeMap.get(o2).size()).compareTo(typeMap.get(o1).size()));
		return genes;
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
			}
		}
	}

	protected void printUniqueAA()
	{
		Set<String> sites = new HashSet<>();
		for (String gene : typeMap.keySet())
		{
			sites.addAll(typeMap.get(gene).keySet());
		}
		TermCounter tc = new TermCounter();
		for (String site : sites)
		{
			tc.addTerm(site.substring(0, 1));
		}
		tc.print();
	}

	public void fillInMissingEffect(Collection<ProteomicsFileRow> datas, int proximityThreshold)
	{
		for (ProteomicsFileRow data : datas)
		{
			if (data.effect != null) continue;
			if (data.sites == null || data.sites.isEmpty()) continue;

			Set<Integer> found = getEffects(data, proximityThreshold);
			data.effect = aggregateEffects(found);
		}
	}

	public Set<Integer> getEffects(ProteomicsFileRow data, int proximityThreshold)
	{
		Set<Integer> found = new HashSet<>();

		for (String gene : data.sites.keySet())
		{
			for (String site : data.sites.get(gene))
			{
				Integer e = getEffect(gene, site);
				if (e != null) found.add(e);
			}
		}

		if (found.isEmpty() && proximityThreshold > 0)
		{
			for (String gene : data.sites.keySet())
			{
				for (String site : data.sites.get(gene))
				{
					Integer e = getClosestEffect(gene, site, proximityThreshold);
					if (e != null) found.add(e);
				}
			}
		}

		return found;
	}

	protected ProteomicsFileRow.SiteEffect aggregateEffects(Set<Integer> found)
	{
		if (found.contains(1))
		{
			if (found.contains(-1)) return ProteomicsFileRow.SiteEffect.COMPLEX;
			else return ProteomicsFileRow.SiteEffect.ACTIVATING;
		}
		else if (found.contains(-1))
		{
			return ProteomicsFileRow.SiteEffect.INHIBITING;
		}
		else if (!found.isEmpty()) return ProteomicsFileRow.SiteEffect.COMPLEX;
		return null;
	}
}
