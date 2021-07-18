package org.panda.resource.siteeffect;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;

/**
 * Serves the Signor database.
 *
 * @author Ozgun Babur
 */
public class CustomSiteEffectServer extends SiteEffectServer
{
	public boolean load(String file) throws IOException
	{
		typeMap = new HashMap<>();
		Feature mod = Feature.PHOSPHORYLATION;
		typeMap.put(mod, new HashMap<>());

		Files.lines(Paths.get(file)).map(l -> l.split("\t")).forEach(t ->
		{
			String gene = t[0];
			String site = t[1];
			Integer effect = Integer.valueOf(t[2]);

			if (!typeMap.get(mod).containsKey(gene)) typeMap.get(mod).put(gene, new HashMap<>());
			typeMap.get(mod).get(gene).put(site, effect);
		});

		return true;
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[0];
	}

	@Override
	public boolean load() throws IOException
	{
		return true;
	}
}
