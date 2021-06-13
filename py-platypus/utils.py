import pickle


def save_pickle(name, data):
    '''save some data to a pickle file'''
    file_name = name + ".p" 
    pickle.dump(data, open(file_name, "wb"))
    return
