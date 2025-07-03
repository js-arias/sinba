#include <cpp11.hpp>
#include <cpp11armadillo.hpp>
#include "cpp11/matrix.hpp"
#include "cpp11/doubles.hpp"
#include <limits>
using namespace cpp11;
using namespace arma;
namespace writable = cpp11::writable;

// mult_by_time multiplies a Q matrix by time.
doubles_matrix<> mult_by_time(writable::doubles_matrix<> Q, double t)
{
	for (int i = 0; i < Q.nrow(); i++)
	{
		for (int j = 0; j < Q.ncol(); j++)
		{
			Q(i, j) *= t;
		}
	}
	return Q;
}

// conditional calculates the conditional likelihood
// of a branch segment.
doubles conditional(doubles to, doubles_matrix<> Q, double time)
{
	double mx = to[0];
	for (int i = 1; i < to.size(); i++)
	{
		if (to[i] > mx)
		{
			mx = to[i];
		}
	}

	mat pm = expmat(as_mat(mult_by_time(Q, time)));

	writable::doubles x(to.size());
	for (int i = 0; i < to.size(); i++)
	{
		double sum = 0;
		for (int j = 0; j < to.size(); j++)
		{
			sum += pm(i, j) * exp(to[j] - mx);
		}
		x[i] = log(sum) + mx;
	}

	return x;
}

// sinba_conditionals calculates the likelihood conditionals
// using a vector anc with parent nodes,
// vector desc with descendant nodes,
// vector st with the activate status of the node,
// vector lengths with the the branch lengths of each edge,
// a matrix with the conditionals of the tips
// the ages (in a brach) of each event,
// and the Q matrices.
[[cpp11::register]]
cpp11::doubles_matrix<> sinba_conditionals(integers anc, integers desc,
					   integers st, doubles lengths,
					   writable::doubles_matrix<> cond,
					   double first_age, double second_age,
					   doubles_matrix<> root_Q, doubles_matrix<> semi_Q, doubles_matrix<> Q)
{
	// edges are sorted from the root
	for (int i = desc.size() - 1; i >= 0; i--)
	{
		int n = desc[i];
		int a = anc[i];

		writable::doubles to(cond.ncol());
		for (int j = 0; j < cond.ncol(); j++)
		{
			to[j] = cond(n, j);
		}

		doubles from;
		switch (st[a])
		{
		case 0:
			if (st[n] == 2)
			{
				from = conditional(to, semi_Q, lengths[i]);
				break;
			}
			from = conditional(to, root_Q, lengths[i]);
			break;
		case 1:
			from = conditional(to, semi_Q, lengths[i]);
			break;
		default:
			from = conditional(to, Q, lengths[i]);
		}

		for (int j = 0; j < cond.ncol(); j++)
		{
			cond(a, j) += from[j];
		}
	}
	return cond;
}