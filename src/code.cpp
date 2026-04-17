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
			double p = pm(i, j);
			if (p <= 0)
			{
				// in some particular cases
				// there are values smaller than 0
				// (because of rounding)
				continue;
			}
			sum += pm(i, j) * exp(to[j] - mx);
		}
		x[i] = log(sum) + mx;
	}

	return x;
}

// full_conditional calculates the likelihood conditionals
// for a whole tree
// using a vector anc with parent nodes,
// vector desc with descendant nodes,
// vector lengths with the the branch lengths of each edge,
// a matrix with the conditionals of the tips,
// and the Q matrices.
[[cpp11::register]]
cpp11::doubles_matrix<> full_conditionals(integers anc, integers desc,
					  doubles lengths,
					  writable::doubles_matrix<> cond,
					  doubles_matrix<> Q)
{
	// edges are sorted from the root
	for (int i = desc.size() - 1; i >= 0; i--)
	{
		int n = desc[i];
		int a = anc[n];

		writable::doubles to(cond.ncol());
		for (int j = 0; j < cond.ncol(); j++)
		{
			to[j] = cond(n, j);
		}

		doubles from = conditional(to, Q, lengths[n]);
		for (int j = 0; j < cond.ncol(); j++)
		{
			cond(a, j) += from[j];
		}
	}
	return cond;
}

// FitzJohn prior
doubles fitzJohnPrior(doubles to, doubles_matrix<> PI_mat, int size)
{
	double mx = to[0];
	for (int i = 1; i < to.size(); i++)
	{
		if (to[i] > mx)
		{
			mx = to[i];
		}
	}
	double sum = 0;
	for (int i = 0; i < to.size(); i++)
	{
		sum += exp(to[i] - mx);
	}

	writable::doubles pi(size);
	for (int i = 0; i < pi.size(); i++)
	{
		pi[i] = 0;
	}

	for (int i = 0; i < PI_mat.ncol(); i++)
	{
		for (int j = 0; j < PI_mat.ncol(); j++)
		{
			int pos = int(PI_mat(i, j));
			if (pos == 0)
			{
				continue;
			}
			pos--;
			pi[pos] += exp(to[j] - mx) / sum;
		}
	}
	return (pi);
}

// birthEvent updates the conditionals on a birth event.
doubles birth_event(doubles to, doubles_matrix<> PI_mat, doubles pi)
{
	double sum_pi = 0;
	for (int i = 0; i < pi.size(); i++)
	{
		sum_pi += pi[i];
	}
	if (sum_pi < 1e-6)
	{
		pi = fitzJohnPrior(to, PI_mat, pi.size());
	}

	double mx = to[0];
	for (int i = 1; i < to.size(); i++)
	{
		if (to[i] > mx)
		{
			mx = to[i];
		}
	}

	writable::doubles x(to.size());
	for (int i = 0; i < to.size(); i++)
	{
		double sum = 0;
		for (int j = 0; j < to.size(); j++)
		{
			int pos = int(PI_mat(i, j));
			if (pos == 0)
			{
				continue;
			}
			pos--;
			sum += pi[pos] * exp(to[j] - mx);
		}
		x[i] = log(sum) + mx;
	}

	return (x);
}

// sinba_conditionals calculates the likelihood conditionals
// using a vector anc with parent nodes,
// vector desc with descendant nodes,
// vector st with the activate status of the node,
// vector lengths with the the branch lengths of each edge,
// a matrix with the conditionals of the tips,
// the ages (in a brach) of each event,
// and the Q matrices.
[[cpp11::register]]
cpp11::doubles_matrix<> sinba_conditionals(integers anc, integers desc,
					   integers st,
					   doubles lengths,
					   writable::doubles_matrix<> cond,
					   double first_age, double second_age,
					   doubles_matrix<> semi_Q, doubles_matrix<> Q,
					   doubles act_pi, doubles semi_pi,
					   doubles_matrix<> act_PI_mat, doubles_matrix<> semi_PI_mat)
{
	// edges are sorted from the root
	for (int i = desc.size() - 1; i >= 0; i--)
	{
		int n = desc[i];
		int a = anc[n];

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
				// two events in the same branch
				// first we calculate the conditional for the birth of the full process.
				double len = lengths[n] - second_age;
				doubles tmp2 = conditional(to, Q, len);

				// update the conditionals with the event
				tmp2 = birth_event(tmp2, act_PI_mat, act_pi);

				// then we calculate the conditional for the birth of the semi-active process.
				len = second_age - first_age;
				if (len < 1e-6)
				{
					// simultaneous birth
					// update the conditionals with the event
					tmp2 = birth_event(tmp2, semi_PI_mat, semi_pi);

					// we move directly to the inactive state
					from = tmp2;
					break;
				}
				doubles tmp1 = conditional(tmp2, semi_Q, len);

				// update the conditionals with the event
				tmp1 = birth_event(tmp1, semi_PI_mat, semi_pi);

				// finally we calculate the conditional for the inactive process.
				from = tmp1;
				break;
			}
			if (st[n] == 1)
			{
				// calculate the conditional for the birth of the semi-active process.
				double len = lengths[n] - first_age;
				doubles tmp1 = conditional(to, semi_Q, len);

				// update the conditionals with the event
				tmp1 = birth_event(tmp1, semi_PI_mat, semi_pi);

				// then we calculate the conditional for the inactive process.
				from = tmp1;
				break;
			}

			// inactive process.
			from = to;
			break;
		case 1:
			if (st[n] == 2)
			{
				// calculate the conditional for the birth of the full process.
				double len = lengths[n] - second_age;
				doubles tmp2 = conditional(to, Q, len);

				// update the conditionals with the event
				tmp2 = birth_event(tmp2, act_PI_mat, act_pi);

				// then we calculate the conditional for the semi-active process.
				from = conditional(tmp2, semi_Q, second_age);
				break;
			}

			// semi-active process.
			from = conditional(to, semi_Q, lengths[n]);
			break;
		default:
			// full process.
			from = conditional(to, Q, lengths[n]);
		}

		for (int j = 0; j < cond.ncol(); j++)
		{
			cond(a, j) += from[j];
		}
	}
	return cond;
}

// sinba_simultaneous calculates the likelihood conditionals
// using a vector anc with parent nodes,
// vector desc with descendant nodes,
// vector st with the activate status of the node,
// vector lengths with the the branch lengths of each edge,
// a matrix with the conditionals of the tips,
// the ages (in a brach) of each event,
// and the Q matrices.
[[cpp11::register]]
cpp11::doubles_matrix<> sinba_simultaneous(integers anc, integers desc,
					   integers st,
					   doubles lengths,
					   writable::doubles_matrix<> cond,
					   double birth_age,
					   doubles_matrix<> Q,
					   doubles act_pi,
					   doubles_matrix<> act_PI_mat)
{
	// edges are sorted from the root
	for (int i = desc.size() - 1; i >= 0; i--)
	{
		int n = desc[i];
		int a = anc[n];

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
				// events in the same branch
				// first we calculate the conditional for the birth of the full process.
				double len = lengths[n] - birth_age;
				doubles tmp2 = conditional(to, Q, len);

				// update the conditionals with the event
				tmp2 = birth_event(tmp2, act_PI_mat, act_pi);

				// finally we calculate the conditional for the inactive process.
				from = tmp2;
				break;
			}

			// inactive process.
			from = to;
			break;
		default:
			// full process.
			from = conditional(to, Q, lengths[n]);
		}

		for (int j = 0; j < cond.ncol(); j++)
		{
			cond(a, j) += from[j];
		}
	}
	return cond;
}

// birth_conditional calculates the normalized conditional likelihood
// of a birth for a lineage of length l,
// given a Q matrix
// and obs states.
[[cpp11::register]]
cpp11::doubles_matrix<> birth_conditional(double l, doubles_matrix<> Q, doubles obs,
					  doubles_matrix<> m_PI)
{
	if (l > 1e-6)
	{
		obs = conditional(obs, Q, l);
	}
	double mx = obs[0];
	for (int i = 1; i < obs.size(); i++)
	{
		if (obs[i] > mx)
		{
			mx = obs[i];
		}
	}

	writable::doubles_matrix<> from(obs.size(), obs.size());
	for (int i = 0; i < from.nrow(); i++)
	{
		for (int j = 0; j < from.ncol(); j++)
		{
			from(i, j) = m_PI(i, j) * exp(obs[j] - mx);
		}
	}
	return from;
}