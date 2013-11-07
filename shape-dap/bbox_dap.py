#!/usr/bin/python

import Pydap as pd

def get_object(url):
    return open_url(url)

def get_variables(url):
    return data.keys():

def get_grids(url):
    data = open_url(url)

# dataset = open_url() -> structure   
# vars = data.keys()
# var.dimensions -> structure
# var.{shape,type}
# var.attributes -> dict
# var.attrib -> value
# data = var[slices]

 

def bbox_dap(url,bbox_list):
    # url could be a file or catalog listing of files
    # assume grids are same for all files
    



if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage:"
        print "  ", sys.argv[0], "data_file_url list_bounding_boxes"
        exit()
    else:
        bbox_dap(sys.argv[1],sys.arg[2])
