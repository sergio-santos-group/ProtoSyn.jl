// Generated from /Users/SergioSantos/Dropbox/dev/ProtoSyn.jl/selection/Selection.g4 by ANTLR 4.7.1
import org.antlr.v4.runtime.Lexer;
import org.antlr.v4.runtime.CharStream;
import org.antlr.v4.runtime.Token;
import org.antlr.v4.runtime.TokenStream;
import org.antlr.v4.runtime.*;
import org.antlr.v4.runtime.atn.*;
import org.antlr.v4.runtime.dfa.DFA;
import org.antlr.v4.runtime.misc.*;

@SuppressWarnings({"all", "warnings", "unchecked", "unused", "cast"})
public class SelectionLexer extends Lexer {
	static { RuntimeMetaData.checkVersion("4.7.1", RuntimeMetaData.VERSION); }

	protected static final DFA[] _decisionToDFA;
	protected static final PredictionContextCache _sharedContextCache =
		new PredictionContextCache();
	public static final int
		T__0=1, T__1=2, T__2=3, T__3=4, T__4=5, T__5=6, T__6=7, T__7=8, T__8=9, 
		CHAIN=10, IDENTIFIER=11, INDEX=12, NUMBER=13, WS=14;
	public static String[] channelNames = {
		"DEFAULT_TOKEN_CHANNEL", "HIDDEN"
	};

	public static String[] modeNames = {
		"DEFAULT_MODE"
	};

	public static final String[] ruleNames = {
		"T__0", "T__1", "T__2", "T__3", "T__4", "T__5", "T__6", "T__7", "T__8", 
		"CHAIN", "INT", "CHAR", "IDENTIFIER", "INDEX", "NUMBER", "WS"
	};

	private static final String[] _LITERAL_NAMES = {
		null, "'and'", "'or'", "'('", "')'", "'not'", "'residue'", "'atom'", "'within'", 
		"'of'", "'chain'"
	};
	private static final String[] _SYMBOLIC_NAMES = {
		null, null, null, null, null, null, null, null, null, null, "CHAIN", "IDENTIFIER", 
		"INDEX", "NUMBER", "WS"
	};
	public static final Vocabulary VOCABULARY = new VocabularyImpl(_LITERAL_NAMES, _SYMBOLIC_NAMES);

	/**
	 * @deprecated Use {@link #VOCABULARY} instead.
	 */
	@Deprecated
	public static final String[] tokenNames;
	static {
		tokenNames = new String[_SYMBOLIC_NAMES.length];
		for (int i = 0; i < tokenNames.length; i++) {
			tokenNames[i] = VOCABULARY.getLiteralName(i);
			if (tokenNames[i] == null) {
				tokenNames[i] = VOCABULARY.getSymbolicName(i);
			}

			if (tokenNames[i] == null) {
				tokenNames[i] = "<INVALID>";
			}
		}
	}

	@Override
	@Deprecated
	public String[] getTokenNames() {
		return tokenNames;
	}

	@Override

	public Vocabulary getVocabulary() {
		return VOCABULARY;
	}


	public SelectionLexer(CharStream input) {
		super(input);
		_interp = new LexerATNSimulator(this,_ATN,_decisionToDFA,_sharedContextCache);
	}

	@Override
	public String getGrammarFileName() { return "Selection.g4"; }

	@Override
	public String[] getRuleNames() { return ruleNames; }

	@Override
	public String getSerializedATN() { return _serializedATN; }

	@Override
	public String[] getChannelNames() { return channelNames; }

	@Override
	public String[] getModeNames() { return modeNames; }

	@Override
	public ATN getATN() { return _ATN; }

	public static final String _serializedATN =
		"\3\u608b\ua72a\u8133\ub9ed\u417c\u3be7\u7786\u5964\2\20\177\b\1\4\2\t"+
		"\2\4\3\t\3\4\4\t\4\4\5\t\5\4\6\t\6\4\7\t\7\4\b\t\b\4\t\t\t\4\n\t\n\4\13"+
		"\t\13\4\f\t\f\4\r\t\r\4\16\t\16\4\17\t\17\4\20\t\20\4\21\t\21\3\2\3\2"+
		"\3\2\3\2\3\3\3\3\3\3\3\4\3\4\3\5\3\5\3\6\3\6\3\6\3\6\3\7\3\7\3\7\3\7\3"+
		"\7\3\7\3\7\3\7\3\b\3\b\3\b\3\b\3\b\3\t\3\t\3\t\3\t\3\t\3\t\3\t\3\n\3\n"+
		"\3\n\3\13\3\13\3\13\3\13\3\13\3\13\3\f\3\f\3\r\3\r\3\16\3\16\3\16\7\16"+
		"W\n\16\f\16\16\16Z\13\16\3\17\6\17]\n\17\r\17\16\17^\3\20\6\20b\n\20\r"+
		"\20\16\20c\3\20\3\20\6\20h\n\20\r\20\16\20i\5\20l\n\20\3\20\3\20\5\20"+
		"p\n\20\3\20\6\20s\n\20\r\20\16\20t\5\20w\n\20\3\21\6\21z\n\21\r\21\16"+
		"\21{\3\21\3\21\2\2\22\3\3\5\4\7\5\t\6\13\7\r\b\17\t\21\n\23\13\25\f\27"+
		"\2\31\2\33\r\35\16\37\17!\20\3\2\7\3\2\62;\5\2C\\aac|\4\2GGgg\4\2--//"+
		"\5\2\13\f\16\17\"\"\2\u0086\2\3\3\2\2\2\2\5\3\2\2\2\2\7\3\2\2\2\2\t\3"+
		"\2\2\2\2\13\3\2\2\2\2\r\3\2\2\2\2\17\3\2\2\2\2\21\3\2\2\2\2\23\3\2\2\2"+
		"\2\25\3\2\2\2\2\33\3\2\2\2\2\35\3\2\2\2\2\37\3\2\2\2\2!\3\2\2\2\3#\3\2"+
		"\2\2\5\'\3\2\2\2\7*\3\2\2\2\t,\3\2\2\2\13.\3\2\2\2\r\62\3\2\2\2\17:\3"+
		"\2\2\2\21?\3\2\2\2\23F\3\2\2\2\25I\3\2\2\2\27O\3\2\2\2\31Q\3\2\2\2\33"+
		"S\3\2\2\2\35\\\3\2\2\2\37a\3\2\2\2!y\3\2\2\2#$\7c\2\2$%\7p\2\2%&\7f\2"+
		"\2&\4\3\2\2\2\'(\7q\2\2()\7t\2\2)\6\3\2\2\2*+\7*\2\2+\b\3\2\2\2,-\7+\2"+
		"\2-\n\3\2\2\2./\7p\2\2/\60\7q\2\2\60\61\7v\2\2\61\f\3\2\2\2\62\63\7t\2"+
		"\2\63\64\7g\2\2\64\65\7u\2\2\65\66\7k\2\2\66\67\7f\2\2\678\7w\2\289\7"+
		"g\2\29\16\3\2\2\2:;\7c\2\2;<\7v\2\2<=\7q\2\2=>\7o\2\2>\20\3\2\2\2?@\7"+
		"y\2\2@A\7k\2\2AB\7v\2\2BC\7j\2\2CD\7k\2\2DE\7p\2\2E\22\3\2\2\2FG\7q\2"+
		"\2GH\7h\2\2H\24\3\2\2\2IJ\7e\2\2JK\7j\2\2KL\7c\2\2LM\7k\2\2MN\7p\2\2N"+
		"\26\3\2\2\2OP\t\2\2\2P\30\3\2\2\2QR\t\3\2\2R\32\3\2\2\2SX\5\31\r\2TW\5"+
		"\31\r\2UW\5\27\f\2VT\3\2\2\2VU\3\2\2\2WZ\3\2\2\2XV\3\2\2\2XY\3\2\2\2Y"+
		"\34\3\2\2\2ZX\3\2\2\2[]\5\27\f\2\\[\3\2\2\2]^\3\2\2\2^\\\3\2\2\2^_\3\2"+
		"\2\2_\36\3\2\2\2`b\5\27\f\2a`\3\2\2\2bc\3\2\2\2ca\3\2\2\2cd\3\2\2\2dk"+
		"\3\2\2\2eg\7\60\2\2fh\5\27\f\2gf\3\2\2\2hi\3\2\2\2ig\3\2\2\2ij\3\2\2\2"+
		"jl\3\2\2\2ke\3\2\2\2kl\3\2\2\2lv\3\2\2\2mo\t\4\2\2np\t\5\2\2on\3\2\2\2"+
		"op\3\2\2\2pr\3\2\2\2qs\5\27\f\2rq\3\2\2\2st\3\2\2\2tr\3\2\2\2tu\3\2\2"+
		"\2uw\3\2\2\2vm\3\2\2\2vw\3\2\2\2w \3\2\2\2xz\t\6\2\2yx\3\2\2\2z{\3\2\2"+
		"\2{y\3\2\2\2{|\3\2\2\2|}\3\2\2\2}~\b\21\2\2~\"\3\2\2\2\r\2VX^cikotv{\3"+
		"\2\3\2";
	public static final ATN _ATN =
		new ATNDeserializer().deserialize(_serializedATN.toCharArray());
	static {
		_decisionToDFA = new DFA[_ATN.getNumberOfDecisions()];
		for (int i = 0; i < _ATN.getNumberOfDecisions(); i++) {
			_decisionToDFA[i] = new DFA(_ATN.getDecisionState(i), i);
		}
	}
}