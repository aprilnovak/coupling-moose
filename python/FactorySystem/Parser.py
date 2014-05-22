#!/usr/bin/python

import os, sys, re

import ParseGetPot, Factory
from MooseObject import MooseObject
from Warehouse import Warehouse


class Parser:
  def __init__(self, factory, warehouse):
    self.factory = factory
    self.warehouse = warehouse
    self.params_parsed = set()
    self.params_ignored = set()

  def parse(self, filename):
    try:
      root = ParseGetPot.readInputFile(filename)
    except:
      print "Parse Error: " + filename
      sys.exit(1)

    self._parseNode(root)

    if len(self.params_ignored):
      print "Warning detected during test specification parsing\n  File: " #+ os.path.join(test_dir, filename)
      print '       Ignored Parameter(s): ', self.params_ignored


  def extractParams(self, params, getpot_node):
    full_name = getpot_node.fullName()

    # Populate all of the parameters of this test node
    # using the GetPotParser.  We'll loop over the parsed node
    # so that we can keep track of ignored parameters as well
    local_parsed = set()
    for key, value in getpot_node.params.iteritems():
      self.params_parsed.add(full_name + '/' + key)
      local_parsed.add(key)
      if key in params:
        if params.type(key) == list:
          params[key] = value.split(' ')
        else:
          if re.match('".*"', value):   # Strip quotes
            params[key] = value[1:-1]
          else:
            # Prevent bool types from being stored as strings.  This can lead to the
            # strange situation where string('False') evaluates to true...
            if params.isValid(key) and (type(params[key]) == type(bool())):
              # We support using the case-insensitive strings {true, false} and the string '0', '1'.
              if (value.lower()=='true') or (value=='1'):
                params[key] = True
              elif (value.lower()=='false') or (value=='0'):
                params[key] = False
              else:
                print "Unrecognized (key,value) pair: (", key, ',', value, ")"
                sys.exit(1)

              # Otherwise, just do normal assignment
            else:
              params[key] = value
      else:
        self.params_ignored.add(key)

    # Make sure that all required parameters are supplied
    required_params_missing = params.required_keys() - local_parsed
    if len(required_params_missing):
      print "Error detected during test specification parsing\n  File: " #+ os.path.join(test_dir, filename)
      print '       Required Missing Parameter(s): ', required_params_missing

  # private:
  def _parseNode(self, node):
    if 'type' in node.params:
      moose_type = node.params['type']

      # Get the valid Params for this type
      params = self.factory.validParams(moose_type)

      # Extract the parameters from the Getpot node
      self.extractParams(params, node)

      # Build the object
      moose_object = self.factory.create(moose_type, node.name, params)

      # Put it in the warehouse
      self.warehouse.addObject(moose_object)

    # Loop over the section names and parse them
    for child in node.children_list:
      self._parseNode(node.children[child])