function varargout = savepath(varargin)
% No-op: keeps this study from rewriting the user's pathdef.m.
warning('study:savepath', 'savepath disabled for this study');
varargout = repmat({0}, 1, nargout);
end
