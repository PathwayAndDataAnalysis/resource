package org.panda.resource.siteeffect;

import org.panda.resource.tcga.ProteomicsFileRow;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;

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
}
