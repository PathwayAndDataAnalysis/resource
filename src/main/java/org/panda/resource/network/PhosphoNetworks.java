package org.panda.resource.network;

import org.panda.resource.FileServer;
import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.graph.Graph;
import org.panda.utility.graph.PhosphoGraph;

import java.io.File;
import java.io.IOException;
import java.util.Scanner;

/**
 * Serves the PhosphoNetworks database.
 *
 * @author Ozgun Babur
 */
public class PhosphoNetworks extends FileServer
{
	static PhosphoNetworks instance;

	static PhosphoGraph graph;

	public static PhosphoNetworks get()
	{
		if (instance == null) instance = new PhosphoNetworks();
		return instance;
	}

	public PhosphoGraph getGraph()
	{
		return graph;
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"PhosphoNetworks.txt"};
	}

	@Override
	public String[] getDistantURLs()
	{
		return new String[]{"http://phosphonetworks.org/download/highResolutionNetwork.csv"};
	}

	@Override
	public boolean load() throws IOException
	{
		graph = new PhosphoGraph("PhosphoNetworks", SignedType.PHOSPHORYLATES.getTag());

		Scanner scanner = new Scanner(new File(locateInBase(getLocalFilenames()[0])));
		String target = null;
		while (scanner.hasNextLine())
		{
			String line = scanner.nextLine();

			if (line.startsWith(">"))
			{
				target = line.substring(1);
			}
			else
			{
				String[] token = line.split("\t");
				String source = token[2];
				String site = token[1];

				if (!graph.hasRelation(source, target))
				{
					graph.putRelation(source, target, "", site);
				}
				else
				{
					graph.addSite(source, target, site);
				}
			}
		}

		return true;
	}

	public static void main(String[] args)
	{
		Graph pc = SignedPC.get().getGraph(SignedType.PHOSPHORYLATES);
		Graph pn = get().getGraph();

		pc.printVennIntersections(pn);
	}
}
