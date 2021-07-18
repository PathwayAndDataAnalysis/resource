package org.panda.resource.siteeffect;

import org.panda.resource.FileServer;
import org.panda.resource.tcga.ProteomicsFileRow;
import org.panda.utility.TermCounter;

import java.util.*;

/**
 * Base class for providing site effects.
 *
 * @author Ozgun Babur
 */
public abstract class SiteEffectServer extends FileServer
{
	Map<Feature, Map<String, Map<String, Integer>>> typeMap;

	public Integer getEffect(String gene, String site, Feature mod)
	{
		if (typeMap.containsKey(mod) && typeMap.get(mod).containsKey(gene))
		{
			return typeMap.get(mod).get(gene).get(site);
		}
		return null;
	}

	public Integer getClosestEffect(String gene, String site, Feature mod, int distanceThreshold)
	{
		Map<String, Map<String, Integer>> map = typeMap.get(mod);
		if (map == null) return null;

		if (map.containsKey(gene))
		{
			int s0 = Integer.parseInt(site.substring(1));

			Integer effect = null;
			int closestDist = Integer.MAX_VALUE;

			for (String ss : map.get(gene).keySet())
			{
				Integer eff = map.get(gene).get(ss);

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
				return Integer.valueOf(o1.substring(1)).compareTo(Integer.valueOf(o2.substring(1)));
			}
			catch (NumberFormatException e)
			{
				return 0;
			}
		});
		return list;
	}

	protected List<String> getGenesWithMostSites(Feature mod)
	{
		if (!typeMap.containsKey(mod)) return Collections.emptyList();

		List<String> genes = new ArrayList<>(typeMap.get(mod).keySet());
		Collections.sort(genes, (o1, o2) -> Integer.valueOf(typeMap.get(mod).get(o2).size()).compareTo(typeMap.get(mod).get(o1).size()));
		return genes;
	}

	protected void printSites(String gene)
	{
		System.out.println("Gene: " + gene);
		for (Feature mod : Feature.values())
		{
			if (typeMap.containsKey(mod) && typeMap.get(mod).containsKey(gene))
			{
				for (String site : sortSites(typeMap.get(mod).get(gene).keySet()))
				{
					Integer sign = typeMap.get(mod).get(gene).get(site);
					System.out.print("\tsite: " + site + "\t" + (sign == 1 ? "activating" : sign == -1 ? "inhibiting" : "complex"));
				}
			}
		}
	}

	protected void printUniqueAA(Feature mod)
	{
		if (typeMap.containsKey(mod))
		{
			Set<String> sites = new HashSet<>();
			for (String gene : typeMap.get(mod).keySet())
			{
				sites.addAll(typeMap.get(mod).get(gene).keySet());
			}
			TermCounter tc = new TermCounter();
			for (String site : sites)
			{
				tc.addTerm(site.substring(0, 1));
			}
			tc.print();
		}
	}

	public void fillInMissingEffect(Collection<ProteomicsFileRow> datas, int proximityThreshold)
	{
		for (ProteomicsFileRow data : datas)
		{
			if (data.effect != null && data.effect != ProteomicsFileRow.SiteEffect.COMPLEX) continue;
			if (data.sites == null || data.sites.isEmpty()) continue;

			Set<Integer> found = getEffects(data, data.mod, proximityThreshold);
			ProteomicsFileRow.SiteEffect e = aggregateEffects(found);
			if (e != null) data.effect = e;
		}
	}

	public Set<Integer> getEffects(ProteomicsFileRow data, Feature mod, int proximityThreshold)
	{
		Set<Integer> found = new HashSet<>();

		for (String gene : data.sites.keySet())
		{
			for (String site : data.sites.get(gene))
			{
				Integer e = getEffect(gene, site, mod);
				if (e != null) found.add(e);
			}
		}

		if (found.isEmpty() && proximityThreshold > 0)
		{
			for (String gene : data.sites.keySet())
			{
				for (String site : data.sites.get(gene))
				{
					Integer e = getClosestEffect(gene, site, mod, proximityThreshold);
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

	public Set<String> getAllGenes(Feature mod)
	{
		return typeMap.containsKey(mod) ? typeMap.get(mod).keySet() : Collections.emptySet();
	}

	public Set<String> getSites(String gene, Feature mod)
	{
		if (typeMap.containsKey(mod))
			return typeMap.get(mod).containsKey(gene) ? typeMap.get(mod).get(gene).keySet() : Collections.emptySet();
		return Collections.emptySet();
	}
}
