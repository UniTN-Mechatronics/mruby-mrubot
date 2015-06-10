module Mrubot
  class NewtonRaphson
    attr_accessor :x0, :max_error, :max_iter
    attr_reader :f, :fd
    def initialize
      @x0     = 0.0
      @max_error = 1e-6
      @max_iter  = 100
      @f         = lambda {|x| 2*x**2 - 15*x + 5}
      @fd        = lambda {|x| 4*x - 15}
    end
    alias :guess :x0
    
    def f=(block)
      raise ArgumentError unless block.lambda?
      @f = block
    end

    def fd=(block)
      raise ArgumentError unless block.lambda?
      fd = @block
    end
  end
end