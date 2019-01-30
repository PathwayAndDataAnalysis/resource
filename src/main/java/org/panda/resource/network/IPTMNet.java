package org.panda.resource.network;

import org.panda.resource.FileServer;
import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.graph.PhosphoGraph;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Scanner;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Serves the iPTMNet database.
 *
 * To update, go to http://annotation.dbi.udel.edu/iptmnet_data/ and use the file "ptm.txt". Rename it to iPTMNet.txt.
 *
 * @author Ozgun Babur
 */
public class IPTMNet extends FileServer
{
	static IPTMNet instance;

	static PhosphoGraph graph;

	public static IPTMNet get()
	{
		if (instance == null) instance = new IPTMNet();
		return instance;
	}

	public PhosphoGraph getGraph()
	{
		return graph;
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"iPTMNet.txt"};
	}

	@Override
	public boolean load() throws IOException
	{
		graph = new PhosphoGraph("iPTMNet", SignedType.PHOSPHORYLATES.getTag());

//		Set<String> set = Files.lines(Paths.get(locateInBase(getLocalFilenames()[0]))).map(l -> l.split("\t")[0]).collect(Collectors.toSet());
//		System.out.println("set = " + set);


		Files.lines(Paths.get(locateInBase(getLocalFilenames()[0]))).map(l -> l.split("\t"))
			.filter(t -> t.length > 7)
			.filter(t -> t[0].equals("PHOSPHORYLATION"))
			.filter(t -> t[4].equals("Homo sapiens (Human)"))
			.filter(t -> !t[7].isEmpty())
			.forEach(t ->
		{
			String source = t[7];
			String target = t[3];
			String site = t[5];

			if (!graph.hasRelation(source, target))
			{
				graph.putRelation(source, target, "", site);
			}
			else
			{
				graph.addSite(source, target, site);
			}
		});

		return true;
	}

	public static void main(String[] args)
	{
		PhosphoGraph pn = get().getGraph();
//		pn.printStats();

		Set<String> sites = pn.getSites("SRC", "ITGB3");
		sites.forEach(System.out::println);

//		boolean contains = pn.getUpstream("ITGB3").contains("SRC");
//		System.out.println("contains = " + contains);
	}
}
