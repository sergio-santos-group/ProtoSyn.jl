from src.SelectionListener import SelectionListener

class Scope(object):
    def __init__(self):
        self.scope = {}
        self.counter = 0
    
    def define(self, key, value=None):
        if key not in self.scope:
            if value is None:
                value = '__selection' + str(self.counter)
                self.counter += 1
            self.scope[key] = value
            # self.counter += 1
        return self.scope[key]
    
    def resolve(self, key):
        return self.scope[key]


class CustomSelectionListener(SelectionListener):
    def __init__(self, *args, **kwargs):
        SelectionListener.__init__(self, *args, **kwargs)
        self.scope = Scope()
        self.queue = []
    
    def enqueue(self, cmd):
        self.queue.append(cmd)
    
    def enterSele(self, ctx):
        # self.scope = {}
        print('entering Sele')
        # print(ctx.expr)

    def exitSele(self, ctx):
        print('exiting Sele')
        for cmd in self.queue:
            print(cmd)
        # print(self.scope)
    
    def enterSelectionByID(self, ctx):
        selector = ctx.children[0].getText()
        ID = ctx.INDEX()
        
        varname = self.scope.define(ctx)
        jl = f'{varname} = selection_by_id(molecule, :{selector}, {ID})'

        self.enqueue(jl)
    
    def enterSelectionByName(self, ctx):
        selector = ctx.children[0].getText()
        NAME = ctx.IDENTIFIER()
        
        varname = self.scope.define(ctx)
        # jl = f'{varname} = {selector}_by_name(molecule, {NAME})'
        jl = f'{varname} = selection_by_name(molecule, :{selector}, {NAME})'

        self.enqueue(jl)
    
    def exitSelectionByDistance(self, ctx):
        sele1 = ctx.children[0]
        sele2 = ctx.children[4]
        cutoff = ctx.children[2]

        var1 = self.scope.resolve(sele1)
        var2 = self.scope.resolve(sele2)
        varname = self.scope.define(ctx)
        jl = f'{varname} = select_by_distance(state, {var1}, {var2}, {cutoff.NUMBER()})'
        self.enqueue(jl)
    
    def exitExprAndOrExpr(self, ctx):
        lexpr, rexpr = ctx.expr()
        op = ctx.children[1].getText()
        
        lvar = self.scope.resolve(lexpr)
        rvar = self.scope.resolve(rexpr)
        varname = self.scope.define(ctx)

        jl = f'{varname} = {lvar} {op} {rvar}'
        self.enqueue(jl)
        # pass

    def exitParenthesizedExpr(self, ctx):
        sele = ctx.expr()
        varname = self.scope.define(ctx, self.scope.resolve(sele))
        # print(varname)

    def exitSelectionExpr(self, ctx):
        sele = ctx.selection()
        varname = self.scope.define(ctx, self.scope.resolve(sele))
        # print('exitSelectionExpr',varname,self.scope.resolve(sele))
    