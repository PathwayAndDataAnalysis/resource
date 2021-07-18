package org.panda.resource.siteeffect;

import org.panda.resource.tcga.ProteomicsFileRow;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class SiteEffectCollective
{
	private List<SiteEffectServer> servers;

	public SiteEffectCollective()
	{
		servers = new ArrayList<>();
		servers.add(PhosphoSitePlus.get());
		servers.add(Signor.get());
	}

	public Integer getEffect(String gene, Collection<String> sites, Feature mod)
	{
		Integer eff = null;
		for (String site : sites)
		{
			eff = combineEffect(eff, getEffect(gene, site, mod));
		}
		return eff;
	}

	private Integer combineEffect(Integer soFarEffect, Integer newEffect)
	{
		if (newEffect == null || newEffect == 0) return soFarEffect;
		else if (soFarEffect == null) return newEffect;
		else if (soFarEffect == 0) return 0;
		else if (soFarEffect * newEffect == 1) return soFarEffect;
		else return 0;
	}

	public Integer getEffect(String gene, String site, Feature mod)
	{
		Integer effect = null;
		for (SiteEffectServer server : servers)
		{
			Integer e = server.getEffect(gene, site, mod);
			if (e != null && e != 0) return e;
			else if (effect == null) effect = e;
		}
		return effect;
	}

	public Integer getClosestEffect(String gene, String site, Feature mod, int distanceThreshold)
	{
		Integer effect = null;
		for (SiteEffectServer server : servers)
		{
			Integer e = server.getClosestEffect(gene, site, mod, distanceThreshold);
			if (e != null && e != 0) return e;
			else if (effect == null) effect = e;
		}
		return effect;
	}

	public void fillInMissingEffect(Collection<ProteomicsFileRow> datas, int proximityThreshold)
	{
		for (SiteEffectServer server : servers)
		{
			server.fillInMissingEffect(datas, proximityThreshold);
		}
	}

	public Integer getEffectFromID(String id)
	{
		String[] t = id.split("-");
		String sym = t[0];
		int i = 1;
		while (i < t.length && !isSite(t[i])) sym += "-" + t[i++];

		Set<String> sites = new HashSet<>();
		for (; i < t.length; i++)
		{
			if (isSite(t[i])) sites.add(t[i]);
			else break;
		}

		if (sites.isEmpty()) return null;

		String lett = sites.iterator().next().substring(0, 1);
		Feature mod = lett.equals("S") || lett.equals("T") || lett.equals("Y") ? Feature.PHOSPHORYLATION : lett.equals("A") ? Feature.ACETYLATION : null;

		if (mod == null) throw new RuntimeException("Uncharacterized amino acid: " + lett);

		return getEffect(sym, sites, mod);
	}

	private boolean isSite(String s)
	{
		if (s.length() > 1)
		{
			if (s.startsWith("S") || s.startsWith("T") || s.startsWith("Y") || s.startsWith("K"))
			{
				for (int i = 1; i < s.length(); i++)
				{
					if (!Character.isDigit(s.charAt(i))) return false;
				}
				return true;
			}
		}
		return false;
	}

	public void clearEffects()
	{
		servers.clear();
	}

	public void loadCustomEffectsFromFile(String file) throws IOException
	{
		CustomSiteEffectServer csec = new CustomSiteEffectServer();
		csec.load(file);
		servers.add(csec);
	}

	public static void main(String[] args) throws IOException
	{
		SiteEffectCollective sec = new SiteEffectCollective();
		System.out.println(sec.getEffect("SMAD3", "T8", Feature.PHOSPHORYLATION));
//		writeToFile();
	}

	public static void writeToFile() throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get("/home/ozgun/Documents/Temp/site-effects.txt"));

		Set<String> genes = new HashSet<>();
		SiteEffectCollective sec = new SiteEffectCollective();
		for (SiteEffectServer server : sec.servers)
		{
			for (Feature mod : Feature.values())
			{
				genes.addAll(server.getAllGenes(mod));
			}
		}

		System.out.println("genes = " + genes.size());

		genes.stream().sorted().forEach(gene ->
		{
			for (Feature mod : Feature.values())
			{
				Set<String> sites = new HashSet<>();

				for (SiteEffectServer server : sec.servers)
				{
					sites.addAll(server.getSites(gene, mod));
				}

				sites.stream().sorted().forEach(site -> FileUtil.writeln(gene + "\t" + site + "\t" + mod + "\t" + sec.getEffect(gene, site, mod), writer));
			}
		});

		writer.close();
	}
}
