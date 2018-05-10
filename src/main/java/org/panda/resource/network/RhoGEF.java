package org.panda.resource.network;

import org.panda.resource.FileServer;
import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.PhosphoGraph;

import java.io.IOException;
import java.util.Collections;

/**
 * Serves the PhosphoNetworks database.
 *
 * @author Ozgun Babur
 */
public class RhoGEF extends FileServer
{
	static RhoGEF instance;

	static DirectedGraph graph;

	public static RhoGEF get()
	{
		if (instance == null) instance = new RhoGEF();
		return instance;
	}

	public DirectedGraph getGraph()
	{
		return graph;
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"rho-gef.txt"};
	}

	@Override
	public boolean load() throws IOException
	{
		graph = new PhosphoGraph("rho-gef", SignedType.ACTIVATES_GTPASE.getTag());
		graph.load(locateInBase(getLocalFilenames()[0]), Collections.singleton("activates-gtpase"));

		return true;
	}

	public static void main(String[] args)
	{
		DirectedGraph pn = get().getGraph();
		System.out.println("pn = " + pn);
	}
}
