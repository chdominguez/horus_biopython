from HorusAPI import Plugin

from Blocks.align import align_block
from Blocks.reres import reset_block

plugin = Plugin()

plugin.addBlock(align_block)
plugin.addBlock(reset_block)
