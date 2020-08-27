import os

sequence_map = {}


def readSequenceFile(filepath):

    try:
        open(filepath)
    except NameError:
        print("No input file")
    else:
        if os.path.exists(filepath):
            with open(filepath, 'rb') as f:
                try:
                    content = f.read()
                    content = content.split('\n')
                    # itr = 0
                    for i in range(len(content)):
                        if (content[i]) and (content[i] is not None) and (content[i] is not '') \
                                and (len(content[i]) > 1000):
                            sequence_map[i] = content[i]
                        else:
                            print("Invalid sequence: ", str(i))
                except OSError:
                    print("Error in reading input sequence file.")

    return sequence_map
