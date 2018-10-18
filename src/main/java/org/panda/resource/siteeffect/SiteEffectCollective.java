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

	public Integer getEffect(String gene, String site)
	{
		Integer effect = null;
		for (SiteEffectServer server : servers)
		{
			Integer e = server.getEffect(gene, site);
			if (e != null && e != 0) return e;
			else if (effect == null) effect = e;
		}
		return effect;
	}

	public Integer getClosestEffect(String gene, String site, int distanceThreshold)
	{
		Integer effect = null;
		for (SiteEffectServer server : servers)
		{
			Integer e = server.getClosestEffect(gene, site, distanceThreshold);
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
		writeToFile();
	}

	public static void writeToFile() throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get("/home/ozgun/Documents/Temp/site-effects.txt"));

		Set<String> genes = new HashSet<>();
		SiteEffectCollective sec = new SiteEffectCollective();
		for (SiteEffectServer server : sec.servers)
		{
			genes.addAll(server.getAllGenes());
		}

		System.out.println("genes = " + genes.size());

		genes.stream().sorted().forEach(gene ->
		{
			Set<String> sites = new HashSet<>();

			for (SiteEffectServer server : sec.servers)
			{
				sites.addAll(server.getSites(gene));
			}

			sites.stream().sorted().forEach(site -> FileUtil.writeln(gene + "\t" + site + "\t" + sec.getEffect(gene, site), writer));
		});

		writer.close();
	}
}
