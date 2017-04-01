package org.panda.resource.network;

import org.panda.resource.ResourceDirectory;
import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.Graph;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * Serves the signed and directed version of REACH that contains edges: phosphorylates and dephosphorylates.
 *
 * @author Ozgun Babur
 */
public class SignedPCNoTransfac extends SignedPC
{
	private static SignedPCNoTransfac instance;

	public static SignedPCNoTransfac get()
	{
		if (instance == null) instance = new SignedPCNoTransfac();
		return instance;
	}

	protected String getPrivateDirectory()
	{
		return ResourceDirectory.get() + File.separator + "SignedPCNoTransfac/";
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"curated-signed-false.sif", "SignedPC-woTF.sif", "curated-signed.sif"};
	}

	@Override
	public String[] getDistantURLs()
	{
		return new String[]{GITHUB_REPO_BASE + "curated-signed-false.sif", GITHUB_REPO_BASE + "SignedPC-woTF.sif.gz",
			GITHUB_REPO_BASE + "curated-signed.sif"};
	}

	public static void main(String[] args)
	{
		printNetworkSizes();
		DirectedGraph graph = SignedPCNoTransfac.get().getGraph(SignedType.PHOSPHORYLATES);
		Set<String> set = graph.getDownstream("TP53");
	}

	private static void printNetworkSizes()
	{
		SignedPCNoTransfac spc = new SignedPCNoTransfac();
		for (SignedType type : SignedType.values())
		{
			spc.getGraph(type).printStats();
		}
	}
}
