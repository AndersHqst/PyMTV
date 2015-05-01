from matplotlib.pylab import plot, ylabel, xlabel, savefig, close, title, figtext, grid, yticks
import os

#
#   Scripts to convert the MTV run results
#

def read_run_results(run_results_file):
    """
    Parses the run_result.txt file from MTV
    :param run_results_file:
    :return:
    """
    heuristics = []
    BIC_scores = []
    independent_models = []
    size_of_c = []
    iteration_time = []
    relation_ships = []
    queries = []
    with open(run_results_file) as fd:
        for line in fd:
            if '\n' in line:
                line = line.replace('\n', '')
            chunks = line.split(' ')
            chunks = [c for c in chunks if c != '' and not ('\t' in c)]
            heuristics.append(float(chunks[0]))
            BIC_scores.append(float(chunks[1]))
            queries.append(float(chunks[2]))
            size_of_c.append(int(chunks[3]))
            independent_models.append(int(chunks[4]))
            iteration_time.append(float(chunks[5]))
            relation_ships.append(chunks[6])

    return heuristics[1:], BIC_scores, independent_models[1:], size_of_c[1:], iteration_time[1:], relation_ships[1:], queries[1:]

def plot_BIC_score(BIC_SCORE, path):
    xlabel('|C|')
    ylabel('BIC score')
    grid(True)
    plot(BIC_SCORE)
    savefig(os.path.join(path, 'BIC.png'))
    close()

def plot_heuristic(heuristic, path):
    xlabel('|C|')
    ylabel('h')
    grid(True)
    plot(heuristic)
    savefig(os.path.join(path, 'heuristic.png'))
    close()

def plot_independent_models(independent_models, path):
    xlabel('|C|')
    ylabel('Independent models')
    grid(True)
    yticks([x+1 for x in range(max(independent_models))])
    plot([x+1 for x in range(len(independent_models))],  independent_models)
    savefig(os.path.join(path, 'independent_models.png'))
    close()

def plot_running_time(running_time, path):
    xlabel('|C|')
    ylabel('MTV iteration in secs.')
    grid(True)
    plot([x for x in range(len(running_time))], running_time)
    savefig(os.path.join(path, 'running_time.png'))
    close()

def plot_size_of_c(size_of_c, path):
    xlabel('|C|')
    ylabel('Max model size |Ci|')
    grid(True)
    plot([x+1 for x in range(len(size_of_c))], size_of_c)
    savefig(os.path.join(path, 'size_of_c.png'))
    close()


def plot_run_results(run_result_folder):
    """
    Creates plots of various MTV run results.
    To use this, pass in the folder that was used as output
    folder for MTV with the -o argument.
    :param run_result_folder:
    :return:
    """
    path = os.path.join(run_result_folder, 'run_result.txt')
    heuristics, BIC_scores, independent_models, size_of_c, iteration_time, relation_ships, queries = read_run_results(path)

    plot_BIC_score(BIC_scores, run_result_folder)
    plot_heuristic(heuristics, run_result_folder)
    plot_independent_models(independent_models, run_result_folder)
    plot_size_of_c(size_of_c, run_result_folder)
    plot_running_time(iteration_time, run_result_folder)

