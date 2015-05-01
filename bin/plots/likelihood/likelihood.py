#
# Script to plot likelihood values for a coin toss
#

from matplotlib.pylab import *

def likelihood(outcomes, p_head):
    p_tail = 1 - p_head
    likelihood = 1.0
    for outcome in outcomes:
        if outcome == 0: #tail
            likelihood *= p_tail
        else:
            likelihood *= p_head
    return likelihood

def plot_coin_toss_entropy():
    from scipy.stats import entropy
    ents = []
    xtics_labels = []
    xtic_vals = []
    # plot entropy for probabilities [0;1]
    for p in xrange(0, 101, 1):
        if p % 10 == 0:
            xtic_vals.append(p)
        prob = p / 100.
        if p % 10 == 0:
            xtics_labels.append(prob)
        ents.append(entropy([prob, 1-prob]))
    xlabel('p=Heads')
    ylabel('H(p, 1-p)')
    xticks(xtic_vals, xtics_labels)
    plot(ents)
    show()

plot_coin_toss_entropy()
    
def plotCoinToss(outcomes):
    likelihoods = []

    # plot likelihood for probabilities [0;1]
    for p in xrange(0, 101, 1):
        xtic_vals.append(p)
        p = p / 100.
        xtics_labels.append(p)
        likelihoods.append(likelihood(outcomes, p))
    xlabel('prob. heads')
    ylabel('likelihood')
    plot(likelihoods)
    show()

# toss = [1,1,0,1]
# plotCoinToss(toss)



